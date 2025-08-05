#!/bin/bash

cd /scratch/tsoies-DL/experiments/data/mouse
wget https://www.encodeproject.org/files/ENCFF521IDK/@@download/ENCFF521IDK.fastq.gz &
wget https://www.encodeproject.org/files/ENCFF006WNS/@@download/ENCFF006WNS.fastq.gz &

cd /scratch/tsoies-DL/experiments/data/human
wget https://www.encodeproject.org/files/ENCFF652IAG/@@download/ENCFF652IAG.fastq.gz &
wget https://www.encodeproject.org/files/ENCFF300NKO/@@download/ENCFF300NKO.fastq.gz &
