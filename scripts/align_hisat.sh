#!/bin/bash

T=90

HISAT_DB=$1
GENOME=$2
FASTQ=$3
BAM=$4

# building the hisat2 index
srun -c ${T} --mem 200G -p amd_256M --time 2-0 hisat2-build ${GENOME} ${HISAT_DB}

#align
srun -c ${T} --mem 1000G -p amd_1Tb --time 3-0 hisat2 -p ${T} --rna-strandness RF -q -x ${HISAT_DB} -U ${FASTQ} | samtools view -b -h -Su | samtools sort > ${BAM}
