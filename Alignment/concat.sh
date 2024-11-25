#!/bin/bash
#$ -S /bin/bash
#$ -cwd
cat `find  -name "$1*R1*"|sort`>$1\_1.fastq.gz
cat `find  -name "$1*R2*"|sort`>$1\_2.fastq.gz

