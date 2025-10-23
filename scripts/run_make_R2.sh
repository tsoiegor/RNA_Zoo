#!/bin/bash

srun -c 5 --mem 1800G -p amd_2Tb --time 30-0 /scratch/tsoies-DL/RNA_Zoo/scripts/make_R2.bin /datasets/tsoies-DL/chaika_data/single_cell_1A_cDNA/E250088063_L01_MERGED.fq \
/scratch/tsoies-DL/final_fastq/Carassius_gibelio.fastq \
/scratch/tsoies-DL/final_fastq/Dicrostonyx_torquatus.fastq \
/scratch/tsoies-DL/final_fastq/Ellobius_talpinus.fastq \
/scratch/tsoies-DL/final_fastq/Larus_michahellis.fastq \
/scratch/tsoies-DL/final_fastq/Myotis_brandtii.fastq \
/scratch/tsoies-DL/final_fastq/Panthera_tigris.fastq \
/scratch/tsoies-DL/final_fastq/Sander_lucioperca.fastq \
/scratch/tsoies-DL/final_fastq/Ursus_maritimus.fastq
