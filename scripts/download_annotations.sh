#!/bin/bash

OUT_DIR=$1

#annotaion is not available for:
#Acomys_dimidiatus
#Dicrostonyx torquatus
#Ellobius_talpinus

#Data lost for Perca fluviatilis

#Carassius_gibelio
wget -O ${OUT_DIR}/Carassius_gibelio.gtf.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/023/724/105/GCF_023724105.1_carGib1.2-hapl.c/GCF_023724105.1_carGib1.2-hapl.c_genomic.gtf.gz

#Larus_michahellis
wget -O ${OUT_DIR}/Larus_michahellis.gtf.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/964/199/755/GCF_964199755.1_bLarMic1.1/GCF_964199755.1_bLarMic1.1_genomic.gtf.gz

#Myotis_brandtii
wget -O ${OUT_DIR}/Myotis_brandtii.gtf.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/412/655/GCF_000412655.1_ASM41265v1/GCF_000412655.1_ASM41265v1_genomic.gtf.gz

#Panthera_tigris
wget -O ${OUT_DIR}/Panthera_tigris.gtf.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/018/350/195/GCF_018350195.1_P.tigris_Pti1_mat1.1/GCF_018350195.1_P.tigris_Pti1_mat1.1_genomic.gtf.gz

#Sander_lucioperca
wget -O ${OUT_DIR}/Sander_lucioperca.gtf.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/008/315/115/GCF_008315115.2_SLUC_FBN_1.2/GCF_008315115.2_SLUC_FBN_1.2_genomic.gtf.gz

#Ursus_maritimus
wget -O ${OUT_DIR}/Ursus_maritimus.gtf.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/017/311/325/GCF_017311325.1_ASM1731132v1/GCF_017311325.1_ASM1731132v1_genomic.gtf.gz
