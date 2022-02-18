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
cpu=$3

# Load correct modules
module load bowtie2/${bowtie2_version}
module load samtools/${samtools_version}

# Get correct files
get_samples ${wd} ${fastq_file} ${data_folder}

full_path_fastq="${fastq_files[$((SLURM_ARRAY_TASK_ID-1))]}"
sample_id="${sample_ids[$((SLURM_ARRAY_TASK_ID-1))]}"
fastq=$(basename ${full_path_fastq})

# Check whether script needs to run
if [[ -s "${wd}/processed/${sample_id}/bowtie2/${sample_id}_contaminants" ]]; then
  echo "File already present"
  exit 0
fi

# Create output dirs
cd "${wd}/processed"
mkdir -p "${sample_id}/bowtie2/"

# Run bowtie to remove contaminants
bowtie2 --seedlen=25 \
  --threads ${cpu} \
  --time \
  --un-gz "${wd}/processed/${sample_id}/bowtie2/${sample_id}_filtered.fastq.gz" \
  -x ${bowtie2_index} \
  -U "${wd}/processed/${sample_id}/trimgalore/${fastq}" \
  -S "${wd}/processed/${sample_id}/bowtie2/${sample_id}_contaminants"

# Create contaminant QC file

# Get total number of reads
tot_reads=$(zcat "${wd}/processed/${sample_id}/trimgalore/${fastq}" | echo $((`wc -l`/4)))
echo -e "RiboseQC run for ${sample_id} on `date` \n" >> "$wd/processed/${sample_id}/${sample_id}contaminants.txt"
# Print headers to file
printf '\t%s\t%s\t%s\n' "READ_TYPE" "READS" "PERCENTAGE" >> "$wd/processed/${sample_id}/${sample_id}contaminants.txt"
# Print total no. of reads
printf '%s\t%s\t%s\t%s\n' $sample_id "Total" $tot_reads "100" >> "$wd/processed/${sample_id}/${sample_id}contaminants.txt"

# For each contaminant type, print absolute and relative number of reads
for contaminant_type in tRNA snRNA snoRNA mtDNA rRNA ; do  
  contaminant_reads=`samtools view "${wd}/processed/${sample_id}/bowtie2/${sample_id}_contaminants" | grep -o "$contaminant_type" | wc -l`
  contaminant_percentage=`awk -vn=248 "BEGIN{print(${contaminant_reads}/${tot_reads}*100)}"`
  printf '%s\t%s\t%s\t%s\n' "${sample_id}" "${contaminant_type}" "${contaminant_reads}" "${contaminant_percentage}" >> "$wd/processed/${sample_id}/${sample_id}contaminants.txt"
done

# Count reads that passed filtering
filtered_reads=$(zcat "${wd}/processed/${sample_id}/bowtie2/${sample_id}_filtered.fastq.gz" | echo $((`wc -l`/4)))  
filtered_percentage=`awk -vn=248 "BEGIN{print(${filtered_reads}/${tot_reads}*100)}"`
printf '%s\t%s\t%s\t%s\n\n' "${sample_id}" "Passed" "${filtered_reads}" "${filtered_percentage}" >> "$wd/processed/${sample_id}/${sample_id}contaminants.txt"