#!/bin/bash


dir=$1

for genome in ${dir}/*; do
new_name=$(zcat ${genome} | head -n 1 | awk -v OFS='_' '{print $2,$3}').fna.gz
echo $(basename ${genome}) '-->' ${new_name}
out_dir=$(dirname ${genome})
echo $out_dir
mv ${genome} ${out_dir}/${new_name}
done
