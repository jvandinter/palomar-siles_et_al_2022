#!/bin/bash

######################################################################
#
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Date:    13-12-2021
# Update   13-12-2021
#
######################################################################

set -uo pipefail

# Load parameters from main script
source $1
source $2
gtf=$3
cpu=$4

# Load correct modules
module load subread/${subread_version}

# Create output dir
mkdir -p "${wd}/data/processed/featurecounts"

echo "`date` running FeatureCounts for ${merged_gtf_basename} samples"

# Check whether script needs to run
if [[ -f "${wd}/data/processed/featurecounts/${merged_gtf_basename}.counts" ]]; then
  echo "`date` Count file already present"
fi

# Get correct files
get_samples $wd $fastq_file $data_folder

for i in ${!sample_ids[@]}; do
  bam_files[i]="${wd}/data/processed/star/${sample_ids[i]}/${sample_ids[i]}.Aligned.sortedByCoord.out.bam"
done

echo "`date` Running featureCounts"
featureCounts -v

featureCounts \
  -s 2 \
  -T ${cpu} \
  -t "CDS" \
  -g "gene_id" \
  -J \
  -G ${reference_genome} \
  -a "${gtf}" \
  -o "${wd}/data/processed/featurecounts/${merged_gtf_basename}.counts" \
  ${bam_files[@]}

echo "`date` Finished ${merged_gtf_basename} samples"