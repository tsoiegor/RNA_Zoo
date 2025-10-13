#!/bin/bash

##### SLURM PARAMS #####
CPU=100
MEM=1000G
P="amd_1Tb"
TIME="30-0"

##### RUN PARAMS #####
exec_path=/scratch/tsoies-DL/chaika_scripts/dnbc4tools2.1.3
fastqDIR=/scratch/tsoies-DL/final_fastq
TMP=/scratch/tsoies-DL/tmp
export TMPDIR=${TMP}
genomeDIR=/datasets/tsoies-DL/genomes
oligoDIR=/datasets/tsoies-DL/chaika_data/single_cell_1A_ologo
OUT_DIR=/scratch/tsoies-DL/Output_DNBelab


for genome_annotation in /scratch/tsoies-DL/annotations/*gtf; do

sp_name=$(basename -s .gtf ${genome_annotation})
echo "working with: ${sp_name}"
mkdir ${OUT_DIR}/${sp_name}
genome=${genomeDIR}/${sp_name}.fna.gz
cDNAfastq1=${fastqDIR}/${sp_name}.fastq.R1.fastq
cDNAfastq2=${fastqDIR}/${sp_name}.fastq

### uncompressed genome ###
mkdir ${OUT_DIR}/${sp_name}/genome
uncompressed_genome=${OUT_DIR}/${sp_name}/genome/${sp_name}.fasta
zcat ${genome} > ${uncompressed_genome}

### processing gtf ###
srun -c ${CPU} --mem ${MEM} --time ${TIME} -p ${P} ${exec_path}/dnbc4tools tools mkgtf --ingtf ${genome_annotation} --output ${OUT_DIR}/${sp_name}/${sp_name}_filtered.gtf --type gene_type

### STAR genome ###
srun -c ${CPU} --mem ${MEM} --time ${TIME} -p ${P} ${exec_path}/dnbc4tools tools mkref --genomeDir ${OUT_DIR}/${sp_name}/genome --ingtf ${OUT_DIR}/${sp_name}/${sp_name}_filtered.gtf --fasta ${uncompressed_genome} --threads ${T} --species ${sp_name}

### running analysis ###
srun -c ${CPU} --mem ${MEM} --time ${TIME} -p ${P} ${exec_path}/dnbc4tools rna run \
    --cDNAfastq1 ${cDNAfastq1} \
    --cDNAfastq2 ${cDNAfastq2} \
    --oligofastq1 $(ls -m ${oligoDIR}/*1.fq.gz | sed -e 's/, /,/g') \
    --oligofastq2 $(ls -m ${oligoDIR}/*2.fq.gz | sed -e 's/, /,/g') \
    --genomeDir  ${OUT_DIR}/${sp_name}/genome \
    --name ${sp_name} --threads ${T}
done
