#!/bin/bash

function check_annotation() {

    # Check whether to use a custom GTF / Rannot or use the prebuild
    # reference GTF / Rannot

    reference_annotation=$1
    reference_gtf=$2
    custom_annotation=$3
    custom_gtf=$4

    # Check to see whether custom GTF / Rannot was provided
    if [[ -f ${custom_gtf} ]]; then
        gtf=${custom_gtf}
    else
        gtf=${reference_gtf}
    fi

    if [[ -f ${custom_annotation} ]]; then
        rannot=${custom_annotation}
    else
        rannot=${reference_annotation}
    fi
        annot_name=$(basename ${rannot%.gtf_Rannot})
        annot_name=${annot_name#Homo_sapiens*.*}
}

function get_samples() {

    # Creates an array that contains all samples that will be analysed, including
    # additional arrays that contain the barcode and sample_ID

    wd=$1
    fastq_file=$2
    data_folder=$3

# Either get fastq files from symlinks in ./raw/ or a list from the config file
    fastq_files=()
    if [[ $(ls "${wd}/raw/" | wc -l) -eq 0 ]]; then

        mapfile -t fastq_names < ${fastq_file}

        if [[ "${fastq_names[0]}" =~ "/hpc/pmc_vanheesch/" ]]; then

            for i in ${!fastq_names[@]}; do
            fastq_files[i]=${fastq_names[i]}
            done

        else

            for i in ${!fastq_names[@]}; do
            fastq_files[i]=$(find "${data_folder}" -name ${fastq_names[i]})
            done
        fi

    else

        fastq_names=($(find "${wd}/raw" -maxdepth 1 -name "*.fastq.gz" -exec basename {} \; | sort -u ))

        for i in ${!fastq_names[@]}; do

            fastq_files[i]=$(find "${data_folder}" -name ${fastq_names[i]})

        done

    fi

# Initiate arrays
    sample_ids=()
    barcode_ids=()
    samples=()

# Get sample IDs and barcodes from fastq files
    for i in ${!fastq_files[@]}; do

        sample_ids[i]=$(basename ${fastq_files[i]} | cut -f 1 -d "_")
        barcode_ids[i]=$(basename ${fastq_files[i]} | cut -f 2 -d "_")
        samples[i]=$(basename ${fastq_files[i]})

    done
}