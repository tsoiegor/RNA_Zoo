#!/bin/bash

T=95
GENOME=$1
RM_out=$2
species=$3
TMP=/scratch/tsoies-DL/tmp
mkdir ${RM_out}

echo "decompressing genome..."

zcat ${GENOME} > ${TMP}/$(basename -s .fna.gz ${GENOME}).fasta
GENOME=${TMP}/$(basename -s .fna.gz ${GENOME}).fasta
echo "decompressing finished"

srun -c ${T} --mem 1000G -p amd_1Tb --time 30-0 RepeatMasker -dir ${RM_out} -pa ${T} -species ${species} -xsmall ${GENOME}
