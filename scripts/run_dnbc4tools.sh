#!/bin/bash

##### SLURM PARAMS #####
CPU=90
MEM=1000G
P="amd_2Tb"
TIME="30-0"

##### RUN PARAMS #####
exec_path=/scratch/tsoies-DL/dnbc4tools2.1.3
fastqDIR=/scratch/tsoies-DL/final_fastq
TMP=/scratch/tsoies-DL/tmp
genomeDIR=/datasets/tsoies-DL/genomes
oligoDIR=/datasets/tsoies-DL/chaika_data/single_cell_1A_ologo
OUT_DIR=/scratch/tsoies-DL/Output_DNBelab

##### CONFIGURING ENVIRONMENT #####
#export PATH=$PATH:${exec_path}/external/conda/bin:${exec_path}/lib/python/dnbc4tools:${exec_path}/lib/python/dnbc4tools/rna
export TMPDIR=${TMP}
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${exec_path}/external/conda/lib

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
srun -c ${CPU} --mem ${MEM} --time ${TIME} -p ${P} ${exec_path}/dnbc4tools tools mkgtf --ingtf ${genome_annotation} --output ${OUT_DIR}/${sp_name}/${sp_name}_filtered.gtf

### STAR genome ###
srun -c ${CPU} --mem ${MEM} --time ${TIME} -p ${P} ${exec_path}/dnbc4tools rna mkref --genomeDir ${OUT_DIR}/${sp_name}/genome --ingtf ${OUT_DIR}/${sp_name}/${sp_name}_filtered.gtf --fasta ${uncompressed_genome} --threads ${CPU} --species ${sp_name}

#srun -c ${CPU} --mem ${MEM} --time ${TIME} -p ${P} STAR --runThreadN ${CPU} --runMode genomeGenerate --genomeDir ${OUT_DIR}/${sp_name}/genome --genomeFastaFiles ${uncompressed_genome} --sjdbGTFfile ${genome_annotation}


### running analysis ###
srun -c ${CPU} --mem ${MEM} --time ${TIME} -p ${P} ${exec_path}/dnbc4tools rna run \
    --cDNAfastq1 ${cDNAfastq1} \
    --cDNAfastq2 ${cDNAfastq2} \
    --oligofastq1 $(ls -m ${oligoDIR}/*1.fq.gz | sed -e 's/, /,/g' | tr -d '\n') \
    --oligofastq2 $(ls -m ${oligoDIR}/*2.fq.gz | sed -e 's/, /,/g' | tr -d '\n') \
    --genomeDir  ${OUT_DIR}/${sp_name}/genome \
    --name ${sp_name} --threads ${CPU} \
    --outdir ${OUT_DIR}/${sp_name}
done
