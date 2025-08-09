#!/bin/bash
echo "usage: merge_genomes_into_one.sh genome_dir out_dir"
genomes_dir=$1
out_dir=$2

read_with='cat'

echo "reading genomes from ${genomes_dir}"
echo "output will be saved into ${out_dir}"

for genome in ${genomes_dir}/*; do
if [[ ${genome} =~ "gz" ]]; then
read_with='zcat'
fi
done
echo "genome files will be red with ${read_with}"
eval "${read_with} ${genomes_dir}/* > ${out_dir}/merged_genome.fa"



