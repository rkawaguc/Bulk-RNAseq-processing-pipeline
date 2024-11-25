##########################################
##############    PCA
##########################################

### load packages
library ("DESeq2")
library(ggplot2)

########## generate object for plot
dds <- DESeqDataSetFromMatrix (mycounts, colData=targets, design = ~ Group)
vsdd <- vst( dds )

theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"), aspect.ratio=1)


###### here we select top 500, top 1000, and top 2000 genes for the PCA plot

pdf("PCAplot.pdf", width=8,height=8)

pcaData <- plotPCA(vsdd, intgroup=c("Group"), ntop = 500, returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"),digits=2)
ggplot(pcaData, aes(PC1, PC2, colour=group)) +
  geom_point(size=2) + geom_text(label=colnames(vsdd),size=2, hjust=0, vjust=0.5)+
  ggtitle("top 500 most variable genes") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + theme+
  coord_fixed()

pcaData <- plotPCA(vsdd, intgroup=c("Group"), ntop = 1000, returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"),digits=2)
ggplot(pcaData, aes(PC1, PC2, colour=group)) +
  geom_point(size=2) + geom_text(label=colnames(vsdd),size=2, hjust=0, vjust=0.5)+
  ggtitle("top 1000 most variable genes") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + theme+
  coord_fixed()

pcaData <- plotPCA(vsdd, intgroup=c("Group"), ntop = 2000, returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"),digits=2)
ggplot(pcaData, aes(PC1, PC2, colour=group)) +
  geom_point(size=2) + geom_text(label=colnames(vsdd),size=2, hjust=0, vjust=0.5)+
  ggtitle("top 2000 most variable genes") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + theme+ 
coord_fixed()

dev.off()

