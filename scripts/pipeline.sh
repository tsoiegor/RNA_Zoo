#!/bin/bash


merge_genomes () {
    GENOMES_DIR=$1
    MERGED_GENOME_DIR=$2
}

find_fastq_file () {
    FASTQ_DIR=$1
}

align_full_database () {
    PATH_TO_FASTQ=$1
    PATH_TO_GENOME=$2
    OUT_BAM_PATH=$3
}

make_barcode_files () {
    OUT_DIR_WITH_BARCODES=$1

}

split_fastq_by_species () {
    OUT_DIR=$1
    DIR_WITH_BARCODES=$2
}




