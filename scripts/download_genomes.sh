#!/bin/bash

out_dir=$1

#Larus
wget -P ${out_dir} https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/964/199/755/GCF_964199755.1_bLarMic1.1/GCF_964199755.1_bLarMic1.1_genomic.fna.gz &

#Dicrostonyx torquatus
wget -P ${out_dir} https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/028/555/095/GCA_028555095.1_ASM2855509v1/GCA_028555095.1_ASM2855509v1_genomic.fna.gz &

#Ellobius talpinus
wget -P ${out_dir} https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/685/095/GCA_001685095.1_ETalpinus_0.1/GCA_001685095.1_ETalpinus_0.1_genomic.fna.gz &

#Acomys dimidiatus
wget -P ${out_dir} https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/907/164/435/GCA_907164435.1_mAcoDim1_REL_1905/GCA_907164435.1_mAcoDim1_REL_1905_genomic.fna.gz &

#Ursus maritimus
wget -P ${out_dir} https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/017/311/325/GCF_017311325.1_ASM1731132v1/GCF_017311325.1_ASM1731132v1_genomic.fna.gz &

#Sander lucioperca
wget -P ${out_dir} https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/008/315/115/GCF_008315115.2_SLUC_FBN_1.2/GCF_008315115.2_SLUC_FBN_1.2_genomic.fna.gz &

#Perca fluviatilis
wget -P ${out_dir} https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/010/015/445/GCF_010015445.1_GENO_Pfluv_1.0/GCF_010015445.1_GENO_Pfluv_1.0_genomic.fna.gz &

#Carassius gibelio
wget -P ${out_dir} https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/023/724/105/GCF_023724105.1_carGib1.2-hapl.c/GCF_023724105.1_carGib1.2-hapl.c_genomic.fna.gz &

#Myotis	brandtii
wget -P ${out_dir} https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/412/655/GCF_000412655.1_ASM41265v1/GCF_000412655.1_ASM41265v1_genomic.fna.gz &

#Panthera tigris
wget -P ${out_dir} https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/018/350/195/GCF_018350195.1_P.tigris_Pti1_mat1.1/GCF_018350195.1_P.tigris_Pti1_mat1.1_genomic.fna.gz 

