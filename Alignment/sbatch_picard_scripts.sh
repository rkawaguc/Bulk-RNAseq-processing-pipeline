#!/bin/bash
#SBATCH -c 6
#SBATCH --mem 50G
#SBATCH --export=ALL

## Riki Kawaguchi modified on 08/08/14 originally from Neelroop Pariksha
## This script submit a single bam or sam file (aligned via STAR program) in the current folder 
## Output is written  in QC folder 
## first qrugment is either usually Aligned.out.sam or Aligned.out.bam
## secpmd argument is either mm9 or mm10 etc.
## third arugment is whether if you want to keep reordered_reads.sorted.bam default=n

TIMESTAMP=$(date +%m%d%y%H%M%S)
                
## Set locations for .jar files
code=picard-tools-1.118

## First, find the Aligned.out.sam file
bfile=$1
genome=$2
keep=${3:-"n"}
workingdir=$PWD
bamdir=$workingdir

case "$genome" in
	"mm10") echo "Will use mm10 genome.";
		refgenome=mm10_genome.fa;
		refFlat=mm10Pred.refFlat;;
	"mm39") echo "Will use mm39 genome.";
		refgenome=mm39_genome.fa;
		refFlat=mm39.refFlat;;
	"hg38") echo "Will use hg38 genome";
		refgenome=hg38_genome.fa;
		refFlat=hg38.refFlat;;
	"mmul") echo "Will use Mmul genome";
                refgenome=mmul_genome.fa;
                refFlat=mmul.refFlat;;
        "mfas") echo "Will use Mfas genome Ensembl refFlat";
                refgenome=mfas_genome.fa;;
                refFlat=mfas.refFlat;;
	"rn6") echo "Will use rn6 genome";
		refgenome=rn6_genome.fa;
		refFlat=rn6.refFlat;;
	*) echo "The genome specified is not supported";;
esac

bfile1=$(basename $bfile)
bfile=${bfile}
echo Current directory is $workingdir
echo Using $bfile as input
if [ -d QC ]; then
	ls ${bamdir}/QC
	if [ `ls ${workingdir}/QC/*log|wc -l` -gt 0 ];then
		echo "Removing old error and out file...."
		rm ${workingdir}/QC/error*
		rm ${workingdir}/QC/out*
  	fi
	if [ -f ${workingdir}/QC/alignment_stats.txt ] && [ -f ${workingdir}/QC/rnaseq_stats.txt ] && [ -f ${workingdir}/QC/gcbias_stats.txt ] && [ -f ${workingdir}/QC/gcbias_summary.txt ] && [ -f ${workingdir}/QC/duplication_stats.txt ] ;then
		echo "Results file already exist. Do you want to erase? [y/n]"
		read answer
		if [ $answer == "y" ];then 
			ls ${workingdir}/QC
			rm ${workingdir}/QC/alignment_stats.txt
			rm ${workingdir}/QC/rnaseq_stats.txt
			rm ${workingdir}/QC/gcbias_stats.txt
			rm ${workingdir}/QC/gcbias_summary.txt
			rm ${workingdir}/QC/duplication_stats.txt
			rm ${workingdir}/QC/duplication_stats.txt
			rm ${workingdir}/QC/insertSizeHist.txt
			echo "removing result files....."
			ls ${workingdir}/QC
		else
			exit
		fi
	fi

	if [ -e ${workingdir}/QC/reordered_reads.sorted.bam ] ;then
		echo "reordered_reads.sorted.bam already exists."
		echo "Do you want to use this file? [y/n]"
		read answer
		case "$answer" in
			"n") rm ${workingdir}/QC/reordered_reads.sorted.bam
		esac
	fi
else 
	echo "Creating QC directory..."
	mkdir ${workingdir}/QC
fi

echo $PATH
echo $bfile 
echo $bamdir 
echo $workingdir 
echo $refgenome 
echo $refFlat 
echo $keep

sbatch -c 6 --mem 50G, -o ${workingdir}/QC/out_${TIMESTAMP}.log -e ${workingdir}/QC/error_${TIMESTAMP}.log ${code}/../RunPicardScripts.sh ${PATH} ${bfile} ${bamdir} ${workingdir} ${refgenome} ${refFlat} ${keep}

exit

