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
rannot=$3
cpu=$4

# Get correct files
get_samples ${wd} ${fastq_file} ${data_folder}
sample_id="${sample_ids[$((SLURM_ARRAY_TASK_ID-1))]}"

# Load software modules
module load R/${r_version}

# Check whether script needs to run
if [[ -f "${wd}/analysis/${sample_id}/RiboseQC/${sample_id}.html" ]]; then
  echo "File already present"
  exit 0
fi

# Create output dirs
cd "${wd}/analysis"
mkdir -p "${sample_id}/RiboseQC/"

# Use RiboSeQC to generate HTML report of the data
Rscript "${scriptdir}/mrp_riboseqc.R" \
  "${wd}/processed/${sample_id}/star/${sample_id}.Aligned.sortedByCoord.out.bam" \
  "${wd}/analysis/${sample_id}/RiboseQC/${sample_id}" \
  "${rannot}" \
  "${annot_name}" \
  "${pandoc_dir}" \
  "${resource_dir}"