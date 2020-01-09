#!/bin/bash
# for i in {1..80}; do echo ${i}_S${i}_L001_R1_001.fastq.gz ${i}_S${i}_L001_R2_001.fastq.gz; done > forward.reverse.ids

module load fastqc/0.11.7

f=forward.reverse.ids
mkdir -p fastq

while read -r line; do
   fastqc $line -o fastq
done < $f

