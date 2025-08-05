GENOMEDIR=/scratch/tsoies-DL/experiments/genomes/star_mix
sample=$1
R1=${sample}/*R1.fastq
R2=${sample}/*R2.fastq
outFileNamePrefix=/scratch/tsoies-DL/experiments/align/mix


export TMPDIR=/scratch/tsoies-DL/tmp
srun -c 20 --mem 500G --time 5-0 -p amd_256M,amd_1Tb STAR --runThreadN 80 --genomeDir $GENOMEDIR --readFilesIn ${R1} ${R2}\
 --outFileNamePrefix ${outFileNamePrefix} --readFilesCommand "cat" --outSAMtype BAM Unsorted --outFilterType BySJout\
 --alignSJoverhangMin 8 --outFilterMultimapNmax 20\
 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999\
 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20\
 --alignIntronMax 1000000 --alignMatesGapMax 1000000\
 --quantMode TranscriptomeSAM --outSAMattributes NH HI AS NM MD
