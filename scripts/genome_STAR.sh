GENOMEDIR=$1
annotation=$2
srun -c 100 --mem 100G -p amd_256M --time 1-0 STAR --runThreadN 4 --runMode genomeGenerate --genomeDir $GENOMEDIR --genomeFastaFiles $GENOMEDIR/*.fasta --sjdbGTFfile ${annotation} --sjdbOverhang 100
