#!/bin/bash

input_fasta=$1
input_fastq=$2
out_bam=$3

srun -c 40 --mem 200G --time 3-0 -p amd_256M bwa mem -t 40 ${input_fasta} ${input_fastq} | srun -c 2 --mem 10G --time 3-0 -p amd_256M samtools view -b -h -Su -@ 2 | srun -c 10 --mem 100G --time 3-0 -p amd_256M samtools sort -@ 10 -T /scratch/tsoies-DL/tmp/tmp > ${out_bam}
