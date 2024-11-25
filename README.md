This site provides standard RNAseq analysis pipeline used in Geschwind Lab. The scripts were written to execute on Orion HPC cluster, however they should work on most of Linux environments.
![Screenshot 2024-11-20 at 1 01 15 PM](https://github.com/user-attachments/assets/8652ca1c-95eb-4ef7-afc6-9e384c584a55)

## 0. Pre-processing fastq files
Concatenating fastq files if multiple runs or lanes, if necessary

create “sample.txt” file which contains sample name

for i in `find ../ -name "*.fastq.gz"`;do echo `basename $i`|cut -f1 -d"_";done|sort|uniq>sample.txt

run the following script in “Data” directory to concatenate the fastq files

for i in `cat sample.txt`;do sbatch /geschwindlabshares/fygaoprj06/GitHub_RNAseqAlignmentPipeline/Script/Alignment/concat.sh $i;done 

Alternatively, you can concatenate files individually, for example,

cat sample1_lane1_R1.fastq.gz sample1_lane2_R1.fastq.gz > sample1_1.fastq.gz

cat sample1_lane1_R2.fastq.gz sample1_lane2_R2.fastq.gz > sample1_2.fastq.gz

## 1. Align the reads to reference genome 

##### Generating genome using STAR 
(For reference, https://github.com/alexdobin/STAR, https://pmc.ncbi.nlm.nih.gov/articles/PMC3530905)

If this is the first time to run STAR, you need to generate genome files with STAR. mm39 genome is available on /geschwindlabshares/fygaoprj06/RNAseq/GenomeMM39/GenomeDir.
If you need other genome for alignment with STAR, use the following command for genome generation.
 
/share/apps/STAR/2.7.5c/bin/Linux_x86_64/STAR --runMode genomeGenerate --genomeDir <genome directory> genome.fa --runThreadN 6  --runThreadN 4 --sjdbFileChrStartEnd sjdbInfo.txt --sjdbOverhang 100

##### Alignment using STAR 

cd /geschwindlabshares/fygaoprj06/RNAseq/TestData

/geschwindlabshares/fygaoprj06/RNAseq/Script/Alignment/rnaseq_STAR.sh /geschwindlabshares/fygaoprj06/RNAseq/GenomeMM39/GenomeDir

There will be individual folders for each sample, and the aligned file names Aligned.out.sam

##### Alignment summary 

/share/apps/R/base/4.0.2/bin/R

setwd("/geschwindlabshares/fygaoprj06/RNAseq/TestData")

source("/geschwindlabshares/fygaoprj06/GitHub_RNAseqPipeline/Script/Alignment/AlignmentSummary.R",echo=T) 

## 2. Convert sam file to bam by samtools and count reads by HTSeq 
(https://github.com/samtools/samtools) (https://htseq.readthedocs.io/en/release_0.11.1/count.html)

##### Reverse stranded library such as Truseq

cd /geschwindlabshares/fygaoprj06/RNAseq/TestData

for i in `find ./ -name Aligned.out.sam`;do pushd `dirname $i`; sbatch /geschwindlabshares/fygaoprj06/GitHub_RNAseqAlignmentPipeline/Script/Alignment/getcountHT_reversestranded.sh /geschwindlabshares/fygaoprj06/RNAseq/GenomeMM39/Mus_musculus.GRCm39.112.gtf `dirname $i | xargs -I {} basename {}`;popd;done


##### Non-stranded library such as SMARTseq v4

cd /geschwindlabshares/fygaoprj06/RNAseq/TestData

for i in `find ./ -name Aligned.out.sam`;do pushd `dirname $i`; sbatch /geschwindlabshares/fygaoprj06/RNAseq/Script/Alignment/getcountHT_nonstranded.sh /geschwindlabshares/fygaoprj06/RNAseq/GenomeMM39/Mus_musculus.GRCm39.112.gtf `dirname $i | xargs -I {} basename {}`;popd;done


## 3. Alignment QC

Install picard (https://broadinstitute.github.io/picard) if necessary. 

##### Running PicardTools 
cd /geschwindlabshares/fygaoprj06/RNAseq/TestData

for i in `find ./ -maxdepth 2 -name Aligned.out.sam`;do pushd `dirname $i`;/geschwindlabshares/fygaoprj06/RNAseq/Script/Alignment/sbatch_picard_scripts.sh Aligned.out.sam mm39 y;popd;done

To compile the results of PicardTools, run the following script.

/geschwindlabshares/fygaoprj06/RNAseq/Script/Alignment/CompileQCOutput.sh

You should see PDF for AlignQC, GCbias, GCsummary, InsertSize, RNAseqc and TranscriptCoverage.
Check the uniformity of transcriptCoverage and InsertSize distribution. StrandBalance should be close to 50% across the board. Note any samples with abnormal values. 

##### Alignment QCs by samtools and bedtools

Raw count barplot provide a quick look of RNA abundance distribution profile. Identify samples with different profile. QC.R contains the script to produce the plot.

![image](https://github.com/user-attachments/assets/0ea41867-745a-47ae-8cac-a893d4997b48)

QC.R
PCA.R

## 4. Normalization and Expression QC

Removing low abundance genes is cruicial for noise reduction. As read alignment is not 100% accurate, low counts genes are mostly noise. Our default is to remove genes if they are not expressed more than 5 counts for 1/3 of samples. However, this is arbitrary and you can use other criteria such as 1PM (counts per million) for at least 3 samples. 

Normalization is important step to remove read depth and read abundance distribution differences. We routinely use the following methods.

1) Trimmed mean of M-values (TMM)
2) Variance stabiliing transformation (VST)
3) Conditional quantile normalization (CQN)

Check the uniformity of expression profile across the samples using boxplot. 

![image](https://github.com/user-attachments/assets/417eb186-4161-46d2-9c52-a0d322ddbb1f)



## 5. Outlier detection
Outlier samples can be a major source of noise. Removing will not only minimizing noise, also improve statistical power and more accurate results from downstream analysis. We identify outliers based on the following criteria. 

1) Different quality detected from alignment rates, RNA complexity, or other QC matrix from picard output.
2) Z-score for normalized counts >2 or < -2
3) Connectivity Z-score > 3 or < -3

![image](https://github.com/user-attachments/assets/27189bf7-91fd-4c73-857a-6efb19722d53)

For 2) and 3), use one of them iteratively until no samples are outside of the criteria. Do not have samples with conditions that hugely affect gene expression (e.g. IP and Input).
The figure above is generated by QC.R

## 6. Regressing out Artifacts

Identifying and regressing sequencing artifacts is crucial to prevent misleading biological conclusions. We use picard output (matrix) and calculate top pincipal components and calculate correlation against top PCs of gene expression matrix.

Example below shows significant correlation of seqPC1-3 with expPC2-3. After regressing them out, these correlation coefficients reduced but not much for expPC1 and group (actually increase for Group aginst expPC2 and 3. It is important to note that this step is done after removing all outlier samples.

![image](https://github.com/user-attachments/assets/5d4204f0-84b1-4162-b006-778d4dc7fee7)

If you have a large number of samples (e.g.>50), we recommend to do this step using bootstrap method. Bootstrap method will select random samples, and calculate beta then repeat this step a number of times (>100). Average beta is calculated and regress it out. This is more conservative method than regressing out beta using the entire samples in case there are some samples that could strongly affect the results. 

Running following R script will produce the correlation plots and regress out sequence artifacts.
Correlation_analysis_and_regress_artifacts.R

#### Removal of Unwanted Variation 

RemoveUnwantedVariations.R

Now the data are ready for downstream analysis such as:

Differential Expression analysis (edgeR, voom, limma etc.) in R

Weighted Correlation Network Analysis with WGCNA package in R


