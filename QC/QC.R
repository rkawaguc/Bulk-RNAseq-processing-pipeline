
## load library
library(MASS)
library(gplots)
library(DESeq2)
library(edgeR)
library(openxlsx)
library(readxl)

#### read in count matrix and get FPKM and logCPM
annot<-read.csv("mm39_annot_ensembl.csv")

mydata<- read.table("metaReadCount.xls",header=T)
rownames(mydata)<-mydata$gene
dim(mydata)

tail(mydata)

mydata.annot<-merge(annot, mydata, by.x="GeneID", by.y="gene")
dim(mydata.annot)
rownames(mydata.annot)<-mydata.annot$GeneID
###the last 5 rows are not the counts

annot<-mydata.annot[,1:6]
mycounts<- mydata.annot[,8:13]

mycounts<-mycounts[,order(colnames(mycounts),decreasing=F)]
colnames(mycounts)

libsize<-colSums(mycounts)
libsize<-as.numeric(libsize)

### get logCPM
dge.obj<- DGEList(counts=mycounts)
logcpm<-cpm(dge.obj,log=TRUE)

## get FPKM
dge.obj$genes$Length <- as.numeric (mydata.annot$Length)
myfpkm<-rpkm(dge.obj)
logfpkm<-log2(myfpkm+0.5)


### read in sample annotation
### note, the order of samples must match the order of the order of columns of mycounts
targets<-read.csv("targets.csv")
labls<- targets$SampleName

>targets
SampleName  Group
Sample1     IP
Sample2     IP
Sample3     tot
Sample4     tot
Sample5     UB
Sample6     UB

### variable stablizing transformed expression
dds <- DESeqDataSetFromMatrix (mycounts, colData=targets, design = ~ Group)
vsdd <- vst( dds )
vsd<-assay(vsdd)

##### color-coding for visualization
group<-as.vector(targets$Group)
group[targets$Group=="IP"]<-"blue"
group[targets$Group=="tot"]<-"magenta"
group[targets$Group=="UB"]<-"darkgreen"

################################################
########## (1) Counts, heatmap -- clustering based on correlation
################################################

matcor1<-cor(mycounts, method="pearson", use="complete")

pdf("SamplePearsonClusteringRawReads.pdf", width=8,height=8)
heatmap.2(as.matrix(matcor1), scale="none", margin=c(7,7),trace="none",col=heat.colors(40), ColSideColors=group, cexCol=0.9,cexRow=0.9, main="Pearson correlation on raw read counts", keysize=1)
dev.off()

###### dendrogram
pdf("dendro_RawRead.pdf",width=10,height=10)
cluster1=hclust(as.dist(1-matcor1),method="average")
plot(cluster1,cex=1.2,labels=labls, xlab="", main="dendrogram, counts")
dev.off()


########## check outliers, by calculating the mean of correlation
IAC1<-matcor1
meanIAC1=apply(IAC1,2,mean)
sdCorr=sd(meanIAC1)
numbersd=(meanIAC1-mean(meanIAC1))/sdCorr

pdf("zscore_Counts.pdf",width=10,height=10)
plot(numbersd,type="n", main=" ", ylab="z-score")
text(numbersd, labels=labls, cex=0.8,col="red")
dev.off()

######### top variable genes
standardize <- function(z) {
rowmed <- apply(z, 1, median)
rowmad <- apply(z, 1, mad)
rv <- sweep(z, 1, rowmed)
rv <- sweep(rv, 1, rowmad, "/")
return(rv)
}

rowMadsT <- apply(mycounts,1,mad)
ordT <- order(rowMadsT,decreasing=TRUE)

top1000T <- head(ordT, n=1000)
T1 <- mycounts[top1000T,]
T1dist <- dist(t(standardize(T1)))
T1dism <- as.matrix(T1dist)

pdf(file="Clustering_Counts_top1000mostVariableGenes.pdf", width=8, height=8)
heatmap.2(T1dism, col=rev(heat.colors(50)), scale="none", trace="none", main="Top 1000 Most Variable Genes - Counts", ColSideColors=group, cexRow=0.9, cexCol=0.9,margin=c(7,7),keysize=1)
dev.off()

################################################
######### (2) Variable Stablizing Transformed Expression
################################################

matcor2 <- cor(vsd, method="pearson", use="complete")
 
pdf("SamplePearsonClustering_VSD.pdf", width=8,height=8) 
heatmap.2(matcor2,scale="none", trace="none",col=heat.colors(40), ColSideColors=group, cexCol=0.9,cexRow=0.9, margin=c(7,7), main=" variance stabilized transformed expression",keysize=1)
dev.off() 

IAC2<-matcor2
pdf("dendro_VSD.pdf",width=10,height=10)
cluster2=hclust(as.dist(1-IAC2),method="average")
plot(cluster2,cex=1.2,labels=labls, xlab="", main="dendrogram, VSD")
dev.off()

meanIAC2=apply(matcor2,2,mean)
sdCorr=sd(meanIAC2)
numbersd=(meanIAC2-mean(meanIAC2))/sdCorr

pdf("zscore_VSD.pdf",width=10,height=10)
plot(numbersd, type="n", ylab="z-score", main="Variance stabilized data z -score")
text(numbersd, labels=colnames(vsd), cex=0.8,col="red")
dev.off()

######### top variable genes
rowMadsT <- apply(vsd,1,mad)
ordT <- order(rowMadsT,decreasing=TRUE)
top1000T <- head(ordT, n=1000)
T2 <- vsd [top1000T,]
T2dist <- dist(t(standardize(T2)))
T2dism <- as.matrix(T2dist)

pdf(file="SampleClustering_VSD_top1000mostVariableGenes.pdf", width=8, height=8)	
heatmap.2(T2dism, col=rev(heat.colors(50)),trace="none", scale="none", main="Top 1000 Most Variable Genes - VSD ", margin=c(7,7), ColSideColors=group, cexRow=0.9, cexCol=0.9,keysize=1)
dev.off()


################################################
######### (3) FPKM
################################################

matcor3 <- cor(logfpkm, method="pearson", use="complete")
  
pdf("SamplePearsonClustering_FPKM.pdf", width=8,height=8)
heatmap.2(matcor3,scale="none", trace="none",col=heat.colors(40), ColSideColors=group, cexCol=0.9,cexRow=0.9, margin=c(7,7), main="Pearson correlation on FPKM",keysize=1)
dev.off()

 
pdf("dendro_FPKM.pdf",width=10,height=10)
IAC3<-matcor3
cluster3=hclust(as.dist(1-IAC3),method="average")
plot(cluster3,cex=1.2,labels=labls, xlab="", main="dendrogram, FPKM")
dev.off()

meanIAC3=apply(matcor3,2,mean)
sdCorr=sd(meanIAC3)
numbersd=(meanIAC3-mean(meanIAC3))/sdCorr

pdf("zscore_FPKM.pdf",width=10,height=10)
plot(numbersd, type="n", ylab="z-score", main="FPKM, z -score")
text(numbersd, labels=colnames(vsd), cex=0.8,col="red")
dev.off()

rowMadsT <- apply(logfpkm,1,mad)
ordT <- order(rowMadsT,decreasing=TRUE)

top1000T <- head(ordT, n=1000)
T3 <- logfpkm [top1000T,]
T3dist <- dist(t(standardize(T3)))
T3dism <- as.matrix(T3dist)

pdf(file="SampleClustering_FPKM_top1000mostVariableGenes.pdf", width=8, height=8)	
heatmap.2(T3dism, col=rev(heat.colors(50)),trace="none", scale="none", main="Top 1000 Most Variable Genes - FPKM ", margin=c(7,7), ColSideColors=group, cexRow=0.9, cexCol=0.9,keysize=1)
dev.off()


################################################
######### (4) logCPM
################################################

matcor4 <- cor(logcpm, method="pearson", use="complete")

pdf("SamplePearsonClustering_logCPM.pdf", width=8,height=8)  
heatmap.2(matcor4,scale="none", trace="none",col=heat.colors(40), ColSideColors=group, cexCol=0.9,cexRow=0.9, margin=c(7,7), main="Pearson correlation on logCPM",keysize=1)
dev.off()


pdf("dendro_logCPM.pdf",width=10,height=10)
IAC4<-matcor4
cluster4=hclust(as.dist(1-IAC4),method="average")
plot(cluster4,cex=1.2,labels=labls, xlab="", main="dendrogram, logCPM")
dev.off()

meanIAC4=apply(matcor4,2,mean)
sdCorr4=sd(meanIAC4)
numbersd4=(meanIAC4-mean(meanIAC4))/sdCorr4

pdf("zscore_logCPM.pdf",width=10,height=10)
plot(numbersd4, type="n", ylab="z-score", main="log CPM, z -score")
text(numbersd4, labels=labls, cex=0.8,col="red")
dev.off()

# variability

rowMadsT <- apply(logcpm,1,mad)
ordT <- order(rowMadsT,decreasing=TRUE)

top1000T <- head(ordT, n=1000)
T4 <- logcpm [top1000T,]
T4dist <- dist(t(standardize(T4)))
T4dism <- as.matrix(T4dist)

pdf(file="SampleClustering_Top_logCPM.pdf", width=8, height=8)	
heatmap.2(T4dism, col=rev(heat.colors(50)),trace="none", scale="none", main="Top 1000 Most Variable Genes - logCPM ", margin=c(7,7), ColSideColors=group,  cexRow=0.9, cexCol=0.9,keysize=1)
dev.off()


################################################
######### (5) Distribution
################################################

######## barplot counts
geneCounts <- data.frame(t(1:6*0))
colnames(geneCounts) <- c("cnts=0", "1<=cnts<10", "10<=cnts<20","20<=cnts<50","50<=cnts<100","100<=cnts<1e6")[6:1]
for(i in 1:ncol(myCounts)){
geneCounts[i,] <- hist(myCounts[,i], breaks=c(0,1,10,20,50,100,10e+100)-0.1, plot=F)$counts[6:1]
}
rownames(geneCounts) <- colnames(mycounts)

pdf(file="rawCounts_Barplot.pdf", width=10, height=10)
barplot(t(geneCounts),axes=FALSE, col=rainbow(6), space=0.01, cex.axis=0.3, las=1, names.arg=1:ncol(myCounts), cex=0.05, xaxt="n",xlab="")
axis(2,at=seq(0,dim(mycounts)[1],5000), las=1, line=-2)  
axis(4,at=seq(0,dim(mycounts)[1],dim(myCounts)[1]/4),c("0","0.25","0.5","0.75", "1.0"), line=-2)
legend("top",c("counts=0","1<=counts<10","10<=counts<20","20<=counts<50","50<=counts<100","100<=counts<1e6")[6:1],fill=rainbow(6))
text(1:length(colnames(mycounts)), par("usr")[3] - 1.5, srt=90, adj=1, labels=colnames(mycounts), cex = 0.7, xpd=TRUE)
dev.off()

######## barplot FPKM
geneCounts <- data.frame(t(1:6*0))
colnames(geneCounts) <- c( "0<=fpkm<1", "1<=fpkm<2","2<=fpkm<5","5<=fpkm<10","10<=fpkm<100","100<=fpkm<1e6")[6:1]
for(i in 1:ncol(vsd3)){
geneCounts[i,] <- hist(myfpkm[,i], breaks=c(0,1,2,5,10,100,10e+100)-0.1, plot=F)$counts[6:1]
}
rownames(geneCounts) <- colnames(myfpkm)
pdf(file="FPKM_Barplot.pdf", width=10, height=10)
barplot(t(geneCounts),axes=FALSE, col=bluered(6), space=0.01, cex.axis=1, las=1, names.arg=1:ncol(vsd2), cex=0.05, xaxt="n",xlab="")
axis(2,at=seq(0, dim(myfpkm)[1], 5000), las=1, line=-2)
axis(4,at=seq(0,dim(myfpkm)[1],dim(myfpkm)[1]/4),c("0","0.25","0.5","0.75", "1.0"), line=-2)
legend("top", c( "0<=fpkm<1", "1<=fpkm<2","2<=fpkm<5","5<=fpkm<10","10<=fpkm<100","100<=fpkm<1e6")[6:1],fill=bluered(6))
text(1:length(colnames(myfpkm)), par("usr")[3] - 1.5, srt=90, adj=1, labels=colnames(myfpkm), cex = 0.7,  xpd=TRUE)
dev.off()

####### box plot
jpeg(file="Boxplot.jpg", width=1200, height=1000)
par(mfrow=c(2,2))
par(cex.axis=0.7)

boxplot( log2(mycounts +0.1) , col=group, main="log2 Counts", xaxt="n",xlab="")
axis(1, labels=FALSE, cex.axis=1)
text(1:length(colnames(mycounts)), par("usr")[3] - 0.25, srt=90, adj=1.2, labels=colnames(mycounts), xpd=TRUE, cex=0.6)

par(cex.axis=0.7)
boxplot(vsd, col= group, main="variance stabilizing transformed expression",xaxt="n",xlab="")
axis(1, labels=FALSE, cex.axis=1)
text(1:length(colnames(vsd)), par("usr")[3] - 0.25, srt=90, adj=1.2, labels=colnames(vsd), xpd=TRUE, cex =0.6)

par(cex.axis=0.7)
boxplot( logfpkm, col= group, main="log2 FPKM", xaxt="n",xlab="")
axis(1, labels=FALSE, cex.axis=1)
text(1:length(colnames(logfpkm)), par("usr")[3] - 0.25, srt=90, adj=1.2, labels=colnames(logfpkm), xpd=TRUE, cex=0.6)

par(cex.axis=0.7)
boxplot( logcpm, col= group, main="log2 CPM", xaxt="n",xlab="")
axis(1, labels=FALSE, cex.axis=1)
text(1:length(colnames(logcpm)), par("usr")[3] - 0.25, srt=90, adj=1.2, labels=colnames(logcpm), xpd=TRUE, cex=0.6)
dev.off()


################################################
######### (6) Multidimensional scaling 
################################################

pdf("MDSplot.pdf", width=10,height=10)
par(mfrow=c(2,2))

plotMDS(vsd, top=1000,labels=labls, pch=19, col=group, cex=0.6, xlab=paste("Dimension 1"), ylab=paste("Dimension 2"),main=" VSD, top 1000 genes")
legend(x="topleft",legend=unique(targets$Group),col=unique(group),fill= unique(group), cex=0.8)

plotMDS(myfpkm, top=1000,labels=labls, pch=19, col=group, cex=0.6, xlab=paste("Dimension 1"), ylab=paste("Dimension 2"),main=" FPKM, top 1000 genes")
legend(x="topleft",legend=unique(targets$Group),col=unique(group),fill= unique(group), cex=0.8)

plotMDS(logcpm, top=1000,labels=labls, pch=19, col=group, cex=0.6, xlab=paste("Dimension 1"), ylab=paste("Dimension 2"),main=" logCPM, top 1000 genes")
legend(x="topleft",legend=unique(targets$Group),col=unique(group),fill= unique(group), cex=0.8)
dev.off()

















