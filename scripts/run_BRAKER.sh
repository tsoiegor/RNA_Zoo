#!/bin/bash

T=48

BAM=$1
GENOME=$2
OUT_DIR=$3
export ORTHODB=/scratch/tsoies-DL/make_annotation/orthoDB/Vertebrata.processed.fa
#export GENEMARK_PATH=/home/tsoies/GeneMark-ETP/bin
#export PROTHINT_PATH=/home/tsoies/ProtHint/bin
export AUGUSTUS_CONFIG_PATH=/scratch/tsoies-DL/make_annotation/BRAKER_OUT/config
#mkdir ${OUT_DIR}
# run BRAKER3
#--prot_seq=${ORTHODB}
srun -c ${T} --mem 500G -p amd_1Tb --time 30-0 singularity exec -B "$(dirname ${AUGUSTUS_CONFIG_PATH}),$(dirname ${OUT_DIR}),$(dirname ${GENOME}),$(dirname ${ORTHODB}),$(dirname ${BAM})" braker3.sif braker.pl --workingdir=${OUT_DIR} \
    --overwrite --genome=${GENOME} --species=$(basename ${OUT_DIR}) \
    --threads ${T} \
    --AUGUSTUS_CONFIG_PATH=${AUGUSTUS_CONFIG_PATH} \
    --bam ${BAM} \
    --prot_seq=${ORTHODB}


#srun -c ${T} --mem 500G -p amd_1Tb --time 30-0 singularity exec -B ${OUT_DIR}:/home/tsoies:rw braker3.sif braker.pl --workingdir=/home/tsoies/sp --overwrite --genome=${GENOME} --species=$(basename ${OUT_DIR}) \
#    --AUGUSTUS_SCRIPTS_PATH=/home/tsoies/miniforge3/pkgs/augustus-3.5.0-pl5321h9716f88_9/bin  \
#    --threads ${T} --AUGUSTUS_BIN_PATH=/beegfs/datasets/home/tsoies/miniforge3/envs/rnaseq/bin \
#    --GENEMARK_PATH=${GENEMARK_PATH} \
#    --bam ${BAM} \
#    --prot_seq=${ORTHODB}
