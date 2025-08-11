#!/bin/bash

echo "usage: genome_to_bowtieIndex.sh FASTA_GENOME OUT_NAME"
CPU=100
MEM_GB=900
FASTA_GENOME=$1
OUT_NAME=$2
echo "SLURM params are: CPU: ${CPU}, MEM: ${MEM_GB}"
echo "genome: ${FASTA_GENOME}, out name: ${OUT_NAME}"

srun -c ${CPU} --mem ${MEM_GB}G --time 1-0 -p amd_1Tb bowtie2-build -f --threads ${CPU} ${FASTA_GENOME} ${OUT_NAME}
