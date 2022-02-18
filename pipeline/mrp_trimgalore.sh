#!/bin/bash

######################################################################
#
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Date: 03-06-2021
#
######################################################################

# Load parameters from main script
source $1
source $2

# Get correct files
get_samples ${wd} ${fastq_file} ${data_folder}

full_path_fastq="${fastq_files[$((SLURM_ARRAY_TASK_ID-1))]}"
sample_id="${sample_ids[$((SLURM_ARRAY_TASK_ID-1))]}"
fastq=$(basename ${full_path_fastq})

# Load correct modules
module load cutadapt/${cutadapt_version}
module load fastqc/${fastqc_version}
module load trimgalore/${trimgalore_version}

# Check whether script needs to run
if [[ -f "${wd}/processed/${sample_id}/trimgalore/${fastq}" ]]; then
  echo "File already present"
  exit 0
fi

# Create output dirs
cd "${wd}/processed"
mkdir -p "${sample_id}/trimgalore/"

# Run trimgalore on both reads
cutadapt --version
fastqc --version

# Change names of trimgalore output
cd "$wd/processed/${sample_id}/trimgalore/"

trim_galore "${full_path_fastq}" \
  --cores 2 \
  --gzip \
  --length 25 \
  --trim-n \
  --fastqc \
  --fastqc_args "--outdir ${wd}/processed/${sample_id}/trimgalore/" \
  --output_dir "$wd/processed/${sample_id}/trimgalore/"

mv "${wd}/processed/${sample_id}/trimgalore/${fastq%.*.*}_trimmed.fq.gz" "${wd}/processed/${sample_id}/trimgalore/${fastq}"

# Calculate trimcounts per paired fastq
tot_reads_1=$(zcat "${full_path_fastq}" | echo $((`wc -l`/4)))

trimmed_reads_1=$(zcat "${wd}/processed/${sample_id}/trimgalore/${fastq}" | echo $((`wc -l`/4)))
trimmed_percentage_1=`awk -vn=248 "BEGIN{print(${trimmed_reads_1}/${tot_reads_1}*100)}"`

# Add read trimming info to run QC file
printf '%s\t%s\t%s\t%s\n' "${sample_id}_1" "Trimmed" $trimmed_reads_1 $trimmed_percentage_1 >> "$wd/processed/${sample_id}/${sample_id}_trim_stats.txt"

