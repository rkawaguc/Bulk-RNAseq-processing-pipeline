######################################
######### Remove Unwanted Variations
######### Adopted from Riki
######################################

###load packages
library(EDASeq)
library(MASS)
library(RUVSeq)
library(DESeq2)

######## read in the count matrix
rawCounts<-read.table("metaReadCount.xls",head=T,row.name=NULL,stringsAsFactors = F)
head(rawCounts);tail(rawCounts)
rownames(rawCounts)<-rawCounts$gene
head(rawCounts)

######### read in sample annotation file
conds<-read.csv("targets.csv",head=T,stringsAsFactors = F)
colnames(conds)=c("sample","group")

rownames(conds)=conds$sample
rawCounts<-rawCounts[,conds$sample]

##### filter lowly expressed genes
temp<-rawCounts
temp$flag<-0
for (i in 1:ncol(rawCounts)){
  temp$flag<-temp$flag+(temp[,i]>=5)
}
table(temp$flag>=3)  
rawCounts<-rawCounts[temp$flag>=3,]

####### generate color coding for visualization
group<-as.factor(conds$group)
cv <- as.matrix(conds)
dim(cv) <- c(1,prod(dim(cv)))
rownames(table(cv))
names <- rownames(table(cv))
col_default <- c("red ","blue","green","orange","bisque4 ","black ","brown ","cyan ","darkgreen ","darkgrey ","darkmagenta ","darkolivegreen ","darkorange ","darkred ","darkslateblue ","darkturquoise ","floralwhite ","greenyellow ","grey ","lightcyan ","lightcyan1 ","lightgreen ","lightsteelblue1 ","lightyellow ","magenta ","mediumpurple3 ","midnightblue ","paleturquoise ","pink ","plum1 ","plum2 ","royalblue ","saddlebrown ","salmon ","sienna3 ","skyblue ","skyblue3 ","steelblue ","tan ","thistle1 ","thistle2 ","turquoise ","violet ","white ","yellowgreen","grey60 ","orangered4 ","brown4 ","darkorange2 ","ivory ")
col <- rainbow(length(names))
names(col)<- names
clab <- matrix(col[as.matrix(conds)],nrow=nrow(conds),ncol=ncol(conds))
colnames(clab) <- colnames(conds) 
rownames(clab)<-rownames(conds)
a=1
for (i in unique(clab[,2])){
  clab[which(clab[,2]==i),2]=col_default[a]
  a=a+1
}
clab


###################################
#######  TMM normalization      
###################################
sizeText=1
design<-model.matrix(~group,data=conds)
dge.obj<-DGEList(counts=rawCounts,group=group)
dge.obj<-calcNormFactors(dge.obj,method="TMM")
dge.obj<-estimateCommonDisp(dge.obj,design)
dge.obj<-estimateGLMCommonDisp(dge.obj,design)
dge.obj<-estimateGLMTagwiseDisp(dge.obj,design)
dge.obj<- estimateGLMTrendedDisp(dge.obj,design)
fit<-glmFit(dge.obj,design)
res<-residuals(fit,type="deviance")
write.csv(dge.obj$pseudo.counts,"TMMnorm.csv")


######################################
pdf("MDSplot_rawCounts.pdf", width=10,height=10)
data <-log2(dge.obj$pseudo.counts+1)
ldat <- dist(t(data))
fit <- isoMDS(ldat, k=2)
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="TMM norm genes with 5 counts at least in 3 samples", type="p",pch=16,cex=1,col=clab[,2],xlim=c(min(x)*1.4,max(x)*1.4))
text(x,y+(max(y)-min(y))/75,labels = row.names(t(data)), cex=sizeText, col=clab[,2])
legend("topright", legend=unique(group), fill=unique(clab[,2]),cex=sizeText)
dev.off()

######################## RUVr k=1
ruvrk1<-RUVr(as.matrix(rawCounts),rownames(rawCounts),k=1,res)
norm.count1<-ruvrk1$normalizedCounts
write.csv(norm.count1,"CountMatrix_RUVrk1.csv")

pdf("MDSplot_RUVrk1.pdf", width=10,height=10)
data <-log2(norm.count1+1)
ldat <- dist(t(data))
fit <- isoMDS(ldat, k=2)
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="RUVr k=1", type="p",pch=16,cex=1,col=clab[,2],xlim=c(min(x)*1.4,max(x)*1.4))
text(x,y+(max(y)-min(y))/75,labels = row.names(t(data)), cex=sizeText, col=clab[,2])
legend("topright", legend=unique(group), fill=unique(clab[,2]),cex=sizeText)
dev.off()


######################## RUVr k=2
ruvrk2<-RUVr(as.matrix(rawCounts),rownames(rawCounts),k=2,res)
norm.count2<-ruvrk2$normalizedCounts
write.csv(norm.count2,"CountMatrix_RUVrk2.csv")

pdf("MDSplot_RUVrk2.pdf", width=10,height=10)
data <-log2(norm.count2+1)
ldat <- dist(t(data))
fit <- isoMDS(ldat, k=2)
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="RUVr k=2", type="p",pch=16,cex=1,col=clab[,2],xlim=c(min(x)*1.4,max(x)*1.4))
text(x,y+(max(y)-min(y))/75,labels = row.names(t(data)), cex=sizeText, col=clab[,2])
legend("topright", legend=unique(group), fill=unique(clab[,2]),cex=sizeText)
dev.off()


######################## RUVr k=3
ruvrk3<-RUVr(as.matrix(rawCounts),rownames(rawCounts),k=3,res)
norm.count3<-ruvrk3$normalizedCounts
write.csv(norm.count3,"CountMatrix_RUVrk3.csv")

pdf("MDSplot_RUVrk3.pdf", width=10,height=10)
data <-log2(norm.count3+1)
ldat <- dist(t(data))
fit <- isoMDS(ldat, k=2)
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="RUVr k=3", type="p",pch=16,cex=1,col=clab[,2],xlim=c(min(x)*1.4,max(x)*1.4))
text(x,y+(max(y)-min(y))/75,labels = row.names(t(data)), cex=sizeText, col=clab[,2])
legend("topright", legend=unique(group), fill=unique(clab[,2]),cex=sizeText)
dev.off()



