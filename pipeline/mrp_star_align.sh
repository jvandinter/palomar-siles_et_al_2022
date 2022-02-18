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
gtf=$3
cpu=$4

# Load correct modules
module load STAR/${star_version}
module load samtools/${samtools_version}

# Get correct files
get_samples ${wd} ${fastq_file} ${data_folder}

sample_id="${sample_ids[$((SLURM_ARRAY_TASK_ID-1))]}"

# Check whether script needs to run
if [[ -s "${wd}/processed/${sample_id}/star/${sample_id}.Aligned.sortedByCoord.out.bam" ]]; then
  echo "File already present"
  exit 0
fi

# Create output dirs
cd "${wd}/processed"
mkdir -p "${sample_id}/star/"

# Map ribo reads
STAR --genomeDir "${star_index_basedir}/29nt" \
  --sjdbGTFfile ${gtf} \
  --runThreadN ${cpu} \
  --runDirPerm All_RWX \
  --twopassMode Basic \
  --readFilesIn \
  "${wd}/processed/${sample_id}/bowtie2/${sample_id}_filtered.fastq.gz" \
  --readFilesCommand zcat \
  --outFilterMismatchNmax 2 \
  --outFilterMultimapNmax 20 \
  --outSAMattributes All \
  --outSAMtype BAM SortedByCoordinate \
  --quantMode GeneCounts \
  --outFileNamePrefix "${wd}/processed/${sample_id}/star/${sample_id}." \
  --limitOutSJcollapsed 10000000 \
  --limitIObufferSize=300000000 \
  --outFilterType BySJout \
  --alignSJoverhangMin 1000 \
  --outTmpKeep None

samtools index -@ ${cpu} "${wd}/processed/${sample_id}/star/${sample_id}.Aligned.sortedByCoord.out.bam"
