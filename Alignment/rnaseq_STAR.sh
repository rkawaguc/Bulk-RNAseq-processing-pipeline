#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V

ref=$1
# There needs to be a "Data:" folder that include the input fastqs

for i in $PWD/Data/*_1.fastq.gz; do nDir=`basename $i _1.fastq.gz`; mkdir $nDir; pushd $nDir; ln -s $i .; end2=`echo $i | sed 's/_1.fastq.gz/_2.fastq.gz/g'`; ln -s $end2 .; sbatch sbatch_STAR.sh *_1.fastq.gz *_2.fastq.gz $ref; popd; done

