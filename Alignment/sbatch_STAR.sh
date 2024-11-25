#!/bin/bash
#SBATCH -c 6
#SBATCH --mem=40G
#SBATCH --export=ALL

end1=$1
end2=$2
ref=$3
echo BEGIN STAR alignment `pwd` `date`
STAR --genomeDir $ref --readFilesIn $end1 $end2 --readFilesCommand zcat --runThreadN 6 --outSAMunmapped Within
echo COMPLETE STAR `pwd` `date`

