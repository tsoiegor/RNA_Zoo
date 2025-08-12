#!/bin/bash

CPU=100
MEM=1000
NODE=amd_1Tb
echo 'USAGE: align_bowtie.sh GENOME FASTQ_LIST OUT_BAM'
GENOME=$1
FASTQ=$2
OUT_BAM=$3
echo "SLURM params CPU: ${CPU};  MEM: ${MEM}GB;  NODE: ${NODE}"
echo "USER ENTERED PARAMS:"
echo " GENOME: ${GENOME}"
echo " FASTQ: ${FASTQ}"
echo " OUT_BAM: ${OUT_BAM}"
echo "alignment in progress.."
export TMPDIR=/scratch/tsoies-DL/tmp

srun -c ${CPU} --mem ${MEM}G -p ${NODE} --time 5-0 bowtie2 --end-to-end --fast --mm -x ${GENOME} -U ${FASTQ} --threads ${CPU} | samtools view -b -h -Su | samtools sort > ${OUT_BAM}

echo "DONE"
