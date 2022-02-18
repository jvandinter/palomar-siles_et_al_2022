#!/bin/bash

#SBATCH -t 72:00:00
#SBATCH --job-name=ribo-sebas

######################################################################
#
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Date: 02-06-2021
#
######################################################################

function usage() {
    cat <<EOF
SYNOPSIS
  mrp_main_script.sh ./path/to/config.config - run orf identifier pipeline
  mrp_main_script.sh help - display this help message
DESCRIPTION
  1. Run TRIMGALORE on ribo-seq reads
  2. Remove contaminants from FASTQ with BOWTIE
  3. Map reads with STAR
  4. Check QC with RIBOSEQC
  5. Run ORFQUANT
  6. Merge BAMs
  7. Run ORFQUANT on merged BAMs

AUTHOR
  Jip van Dinter, MSc
  Damon Hofman, MSc
EOF
}

function info() {
    echo "INFO: $@" >&2
}
function error() {
    echo "ERR:  $@" >&2
}
function fatal() {
    echo "ERR:  $@" >&2
    exit 1
}

# Create a unique prefix for the names for this run of the pipeline. 
# This makes sure that runs can be identified
run=$(uuidgen | tr '-' ' ' | awk '{print $1}')

# Show help message if there was no config file location given on the commandline
if [[ -z $1 ]]; then 

  usage; exit;

fi

# Source all variables from the config file
CONFIG=$1
source ${CONFIG}
source ${scriptdir}/mrp_functions.sh

################################################################################
#
# Find fastq samples in directory
#
################################################################################

get_samples ${wd} ${fastq_file} ${data_folder}

check_annotation ${reference_annotation} ${reference_gtf} ${custom_annotation} ${custom_gtf}

# make sure there are samples
if [[ ${#samples[@]} -eq 0 ]]; then
  fatal "no samples found in ./raw/ or file containing fastq file locations not present"
fi

info "samples: n = ${#samples[@]}"
for sample in ${samples[@]}; do
  info "    $sample"
done

################################################################################
#
# Run the pipeline
#
################################################################################

mkdir -p processed analysis log

echo -e "\n ====== `date` Detect Transcript Isoform Pipeline ====== \n"

echo -e "\n`date` Filtering and trimming ..."
echo -e "====================================================================================== \n"

# 1. TRIMGALORE. Parallel start of all trimgalore jobs to filter for quality
#                with CUTADAPT and output quality reports with FASTQC

trim_jobid=()

trim_jobid+=($(sbatch --parsable \
  --mem=${medium_mem} \
  --cpus-per-task=${high_cpu} \
  --time=24:00:00 \
  --array 1-${#samples[@]}%5 \
  --job-name=${run}.trimgalore \
  --output=log/${run}.trimgalore.%A_%a \
  ${scriptdir}/mrp_trimgalore.sh \
  ${scriptdir}/mrp_functions.sh \
  ${CONFIG}
))

info "trimgalore jobid: ${trim_jobid}"

echo -e "\n`date` Removing contaminants ..."
echo -e "====================================================================================== \n"

# 2. BOWTIE2. Use combination of tRNA, rRNA, snRNA, snoRNA, mtDNA fastas to
#             remove those contaminants from RIBO-seq data. Outputs QC stats
#             to a file per contaminant group.

contaminant_jobid=()

contaminant_jobid+=($(sbatch --parsable \
  --mem=${medium_mem} \
  --cpus-per-task=${medium_cpu} \
  --time=24:00:00 \
  --array 1-${#samples[@]}%5 \
  --job-name=${run}.contaminant \
  --output=log/${run}.contaminant.%A_%a \
  --dependency=aftercorr:${trim_jobid} \
  ${scriptdir}/mrp_remove_contaminants.sh \
  ${scriptdir}/mrp_functions.sh \
  ${CONFIG} \
  ${medium_cpu}
))

info "contaminant jobid: ${contaminant_jobid}"

echo -e "\n`date` Align reads to genome with STAR ..."
echo -e "====================================================================================== \n"

# 3. STAR. Align contaminant-depleted read files to supplied genome and
#          transcriptome. If no new custom transcriptome is supplied, 
#          the normal reference transcriptome is used for 
#          guided assembly.

star_jobid=()

star_jobid+=($(sbatch --parsable \
  --mem=${high_mem} \
  --cpus-per-task=${high_cpu} \
  --time=24:00:00 \
  --array 1-${#samples[@]}%5 \
  --job-name=${run}.star_align \
  --output=log/${run}.star_align.%A_%a \
  --dependency=aftercorr:${contaminant_jobid} \
  ${scriptdir}/mrp_star_align.sh \
  ${scriptdir}/mrp_functions.sh \
  ${CONFIG} \
  ${gtf} \
  ${high_cpu}
))

info "alignment jobid: ${star_jobid}"

echo -e "\n`date` Perform QC with RiboseQC ..."
echo -e "====================================================================================== \n"

# 4. RiboseQC. 

riboseqc_jobid=()

riboseqc_jobid+=($(sbatch --parsable \
  --mem=${medium_mem} \
  --cpus-per-task=${medium_cpu} \
  --time=24:00:00 \
  --array 1-${#samples[@]}%5 \
  --job-name=${run}.riboseqc \
  --output=log/${run}.riboseqc.%A_%a \
  --dependency=aftercorr:${star_jobid} \
  ${scriptdir}/mrp_riboseqc.sh \
  ${scriptdir}/mrp_functions.sh \
  ${CONFIG} \
  ${rannot} \
  ${medium_cpu}
))

info "RiboseQC jobid: ${riboseqc_jobid}"

echo -e "\n`date` Finalise ribo-sebas run ..."
echo -e "====================================================================================== \n"

#10. BASH. Remove redundant files to clear storage after successful run

remove_intermediates_jobid=()

remove_intermediates_jobid+=($(sbatch --parsable \
  --mem=${low_mem} \
  --cpus-per-task=${low_cpu} \
  --time=24:00:00 \
  --job-name=${run}.removal \
  --output=log/${run}.removal \
  --dependency=afterok:${count_occurence_jobid} \
  ${scriptdir}/mrp_cleanup.sh \
  ${CONFIG}
  ))

info "Final jobid: ${remove_intermediates_jobid}"

echo -e "\n ====== `date` Started all jobs! ====== \n"