#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -l h_vmem=4G
#
echo BEGIN `date`
bam=$1
chrHits=`basename $bam .bam`_chrHits.txt
echo -e chrom'\t'`dirname $bam | xargs -I{} basename {}` > $chrHits
echo $bam
/share/apps/samtools/base/1.10/bin/samtools index $bam
/share/apps/samtools/base/1.10/bin/samtools idxstats $bam |head -22|cut -f1,3 | grep -v "_" | awk -F"\t" '{print $1"\t"$2/1000000}' >> $chrHits
echo chromosome count `basename $bam .bam` COMPLETED `date`
