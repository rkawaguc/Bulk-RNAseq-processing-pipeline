#!/bin/bash
#$ -S /bin/bash
#$ -V
#SBATCH -p bigmem
#SBATCH --mem 8G
#$ -cwd

sam=$1
# get random 1M reads from Aligned.out.sam
/share/apps/anaconda3/2020.07/envs/python-2.7/bin/python getRandomSamAllHits.py Aligned.out.sam 1000000

# convert sam to bam
#/share/apps/samtools-1.9/bin/samtools view -Sb sample_any_1000000_Aligned.out.sam > temp.bam
#/share/apps/samtools-1.0/bin/samtools sort temp.bam sorted

/share/apps/samtools/base/1.10/bin/samtools view -Sb sample_any_1000000_Aligned.out.sam > temp.bam
/share/apps/samtools/base/1.10/bin/samtools sort temp.bam> sorted.bam

rm temp.bam


