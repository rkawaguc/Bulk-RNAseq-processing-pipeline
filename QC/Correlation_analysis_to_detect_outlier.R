###################################################
############ Correlation ANalysis
############ Adopted from Riki
###################################################

#### load packages
library(RUVSeq)
library(heatmap.plus)
library(MASS)
library(gplots)
library(DESeq)
library(edgeR)
library(moduleColor)
library(EDASeq)
library(pheatmap)

####### read in stat files from picards
seqQC<-read.table("RNAseqQC.txt",sep=" ",head=T,fill=T,stringsAsFactors = FALSE)
seqQC<-seqQC[,c(1:7,11:21)]
colnames(seqQC)[1]<-"sample"

alignQC<-read.table("RNAseqAlign.txt",sep=" ",head=T,fill=T,stringsAsFactors = FALSE)
alignQC<-alignQC[,c(1,3,4,7,9:12,14:18,21:22)]
colnames(alignQC)[1]="sample"

dupdat<-read.table("RNAseqDuplication.txt",sep=" ",head=TRUE,fill=T,stringsAsFactors = FALSE)
dupdat<-dupdat[,c(1,4,5,7:10)]
colnames(dupdat)[1]="sample"

######### read in count matrix
rawCounts <- read.table("metaReadCount.xls",row.name=1, header=TRUE, fill=TRUE)
rawCounts<-rawCounts[,order(colnames(rawCounts))]

tail(rawCounts)

### remove the last 5 rows that are not counts but alignment stats
if(rownames(tail(rawCounts,1))=="__alignment_not_unique"){
  # remove 5 lines from tail
  rawCounts<-rawCounts[1:(dim(rawCounts)[1]-5),]
}
tail(rawCounts)

####### read in sample annotation file
conds<-read.csv("targets.csv",head=T,row.names=1,stringsAsFactors = F)
colnames(conds)<-c("sample","group")
rownames(conds)<-conds$sample

conds$sample%in%colnames(rawCounts)
rawCounts<-rawCounts[,conds$sample]

seqQC$sample<-gsub("-",".",seqQC$sample)

target<-merge(conds,seqQC,by="sample");dim(conds);dim(seqQC);dim(target)
alignQC$sample<-gsub("-",".",alignQC$sample)
target<-merge(target,alignQC,by="sample");dim(alignQC);dim(target)
dupdat$sample<-gsub("-",".",dupdat$sample)
target<-merge(target,dupdat,by="sample");dim(dupdat);dim(target)

### remove any columns with zero or NA
n=dim(conds)[2]
target<-target[,!(colSums(is.na(target)) > 0)]
temp1<-target[,-c(1:n)]
#temp1<-temp1[,colVars(as.matrix(temp1))!=0]

###scale individual columns after removing outlier
temp1<-t(scale(temp1))
temp1<-as.data.frame(t(temp1))
target<-cbind(target[,1:n],temp1)

write.csv(target,"target2.csv")
rawCounts<-rawCounts[,target$sample]

###############################################
#### Sequencing PC
tar1=target[,c((1+dim(conds)[2]):dim(target)[2])]
thisdat <- t(scale(tar1,scale=F))
PC.metadata <- prcomp(thisdat,center=F);
topPC1 <- PC.metadata$rotation;
varexp <- (PC.metadata$sdev)^2 / sum(PC.metadata$sdev^2)
topvar <- varexp
colnames(topPC1) <- paste("Seq.QC\n",colnames(topPC1)," (",signif(100*topvar[1:6],2),"%)",sep="")
Seq.PC1=as.numeric(topPC1[,1])
Seq.PC2=as.numeric(topPC1[,2])
Seq.PC3=as.numeric(topPC1[,3])
Seq.PC4=as.numeric(topPC1[,4])
Seq.PC5=as.numeric(topPC1[,5])
target$Seq.PC1=Seq.PC1
target$Seq.PC2=Seq.PC2
target$Seq.PC3=Seq.PC3
target$Seq.PC4=Seq.PC4
target$Seq.PC5=Seq.PC5
 

###################################
#          TMM normalization      #
###################################
group<-as.factor(target$group)

temp<-rawCounts
temp$flag<-0
for (i in 1:ncol(rawCounts)){
  temp$flag<-temp$flag+(temp[,i]>=2)
}

table(temp$flag>=2)
rawCounts<-rawCounts[temp$flag>=3,]
filtered<-rawCounts

set<-newSeqExpressionSet(as.matrix(filtered),phenoData=data.frame(group,row.names=colnames(filtered)))
design<-model.matrix(~group,data=pData(set))
dge.obj<-DGEList(counts=filtered,group=xx)
dge.obj<-calcNormFactors(dge.obj,method="TMM")
dge.obj<-estimateCommonDisp(dge.obj,design)
dge.obj<-estimateGLMCommonDisp(dge.obj,design)
dge.obj<-estimateGLMTagwiseDisp(dge.obj,design)
dge.obj<- estimateGLMTrendedDisp(dge.obj,design)
fit<-glmFit(dge.obj,design)
normExpr<-dge.obj$pseudo.counts

###############################
normExpr<-log10(normExpr+1)
datExpr<-normExpr

###generate color-coding forvisualization
cv <- as.matrix(conds)
dim(cv) <- c(1,prod(dim(cv)))
rownames(table(cv))
names <- rownames(table(cv))
col_default <- c("red ","blue","green","orange","bisque4 ","black ","brown ","cyan ","darkgreen ","darkgrey ","darkmagenta ","darkolivegreen ","darkorange ","darkred ","darkslateblue ","darkturquoise ","floralwhite ","greenyellow ","grey ","lightcyan ","lightcyan1 ","lightgreen ","lightsteelblue1 ","lightyellow ","magenta ","mediumpurple3 ","midnightblue ","paleturquoise ","pink ","plum1 ","plum2 ","royalblue ","saddlebrown ","salmon ","sienna3 ","skyblue ","skyblue3 ","steelblue ","tan ","thistle1 ","thistle2 ","turquoise ","violet ","white ","yellowgreen","grey60 ","orangered4 ","brown4 ","darkorange2 ","ivory ")
col <- rainbow(length(names))
names(col)<- names
conds<-conds[target$sample,]
clab <- matrix(col[as.matrix(conds)],nrow=nrow(conds),ncol=ncol(conds))
colnames(clab) <- colnames(conds) 
rownames(clab)<-rownames(conds)
n=dim(conds)[2]
for(k in 1:n){
  a=1;for (i in unique(clab[,k])){
    clab[which(clab[,k]==i),k]=col_default[a]
    a=a+1
  }
}
clab

################# plot MDS before any adjustment
sizeText=0.8
pdf("MDSplot_before_PC_correction.pdf", width=8,height=12)
par(mfrow=c(2,1))
data <- datExpr
ldat <- dist(t(data))
fit <- isoMDS(ldat, k=4)
x <- fit$points[,1]
y <- fit$points[,2]

for(i in 2:dim(conds)[2]){
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="VSD normalized genes with 2 counts for half of samples", type="p",pch=16,cex=1,col=clab[,i],xlim=c(min(x)*1.4,max(x)*1.4))
text(x, y+(max(y)-min(y))/50, labels = row.names(t(data)), cex=sizeText, col=clab[,i])
legend("topleft", legend=unique(as.factor(conds[,i])), fill=unique(clab[,i]),cex=sizeText)
}

x <- fit$points[,3]
y <- fit$points[,4]

for(i in 2:dim(conds)[2]){
  plot(x, y, xlab="Coordinate 3", ylab="Coordinate 4", main="VSD normalized genes with 2 counts for half of samples", type="p",pch=16,cex=1,col=clab[,i],xlim=c(min(x)*1.4,max(x)*1.4))
  text(x, y+(max(y)-min(y))/50, labels = row.names(t(data)), cex=sizeText, col=clab[,i])
  legend("topleft", legend=unique(as.factor(conds[,i])), fill=unique(clab[,i]),cex=sizeText)
}
dev.off()

################# plot MDS before any adjustment for most variable genes
sizeText=0.6
pdf("MDSplot_before_PC_correction_top1000.pdf", width=8,height=12)
rowMadsT <- apply(datExpr,1,mad)
ordT <- order(rowMadsT,decreasing=TRUE)
top1000T <- head(ordT, n=1000)

par(mfrow=c(2,1))
data <- datExpr[top1000T,]
ldat <- dist(t(data))
fit <- isoMDS(ldat, k=4)
x <- fit$points[,1]
y <- fit$points[,2]

for(i in 2:dim(conds)[2]){
  plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="VSD normalized genes with 2 counts for half of samples", type="p",pch=16,cex=1,col=clab[,i],xlim=c(min(x)*1.4,max(x)*1.4))
  text(x, y+(max(y)-min(y))/50, labels = row.names(t(data)), cex=sizeText, col=clab[,i])
  legend("topleft", legend=unique(as.factor(conds[,i])), fill=unique(clab[,i]),cex=sizeText)
}

x <- fit$points[,3]
y <- fit$points[,4]

for(i in 2:dim(conds)[2]){
  plot(x, y, xlab="Coordinate 3", ylab="Coordinate 4", main="VSD normalized genes with 2 counts for half of samples", type="p",pch=16,cex=1,col=clab[,i],xlim=c(min(x)*1.4,max(x)*1.4))
  text(x, y+(max(y)-min(y))/50, labels = row.names(t(data)), cex=sizeText, col=clab[,i])
  legend("topleft", legend=unique(as.factor(conds[,i])), fill=unique(clab[,i]),cex=sizeText)
}
dev.off()

################################################
pdf("PCA_logTMMNorm.pdf",height=20,width=24)

thisdat.HTSC <- t(scale(t(normExpr),scale=F))
PC.HTSC <- prcomp(thisdat.HTSC,center=F);
TopPC1 <- PC.HTSC$rotation[,1:5];
varexp <- (PC.HTSC$sdev)^2 / sum(PC.HTSC$sdev^2)
topvar <- varexp[1:5]
colnames(TopPC1) <- paste("Exp.PC\n",colnames(TopPC1)," (",signif(100*topvar[1:5],2),"%)",sep="")

pairsdat <- data.frame(group=as.numeric(factor(target$group)),
                       #TF=as.numeric(factor(target$TF)),
                       #day=as.numeric(factor(target$day)),
                       #Treatment=as.numeric(factor(target$Treatment)),
                       #Sex=as.numeric(factor(target$Sex)),
                       #Num.sort.x1k=as.numeric(factor(target$Num_Sort.x1k.)),
                       #Batch=as.numeric(factor(target$Batch)),
                       PC1_Sequencing=as.numeric(target$Seq.PC1),
                       PC2_Sequencing=as.numeric(target$Seq.PC2),
                       PC3_Sequencing=as.numeric(target$Seq.PC3),
                       PC4_Sequencing=as.numeric(target$Seq.PC4),
                       PC5_Sequencing=as.numeric(target$Seq.PC5))

cond=labels2colors(as.numeric(factor(target$group)))  ## colors

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) { ## Useful function for comparing multivariate data
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y,use="pairwise.complete.obs",method="pearson"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}


panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) { ## Useful function for comparing multivariate data
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  if (class(x) == "numeric" & class(y) == "numeric") {
    r <- abs(cor(x, y,use="pairwise.complete.obs",method="pearson"))
  } else {
    lmout <- lm(y~x)
    r <- sqrt(summary(lmout)$adj.r.squared)
  }
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
pairs(cbind(TopPC1,pairsdat),col= cond,pch=19,upper.panel = panel.cor,main="Covariates and HTSeqCount Comparison -- |Spearman's rho| correlation values")
dev.off()


################################################
# TopPC1: seqQC, TopPC: Exp PC, TopPC2:mod ExpPC

r2mat_meta = matrix(NA,nrow=5,ncol=dim(conds)[2]-1)
r2mat_seq = matrix(NA,nrow=5,ncol=dim(target)[2]-dim(conds)[2])
rownames(r2mat_meta) <- rownames(r2mat_seq) <- colnames(TopPC1)
datMeta_model = target[,2:n]
datSeq_model=target[,-c(1:n)]

colnames(r2mat_meta) <- colnames(datMeta_model)
colnames(r2mat_seq) <- colnames(datSeq_model)

for(i in c(1:5)){
  for(j in c(1:dim(datMeta_model)[2])){
    to_remove = which(is.na(datMeta_model[,j])==TRUE)
    if(length(to_remove)>0){
      tmp_topPC=TopPC1[-to_remove,i]
      tmp_meta=datMeta_model[-to_remove,j]}
    else{
      tmp_topPC=TopPC1[,i]
      tmp_meta=datMeta_model[,j]  
    }
    mod_mat=model.matrix(~tmp_meta)[,-1]
    mod=summary(lm(tmp_topPC~mod_mat))
    r2mat_meta[i,j]=mod$adj.r.squared
  }
}

for(i in c(1:5)){
  for(j in c(1:dim(datSeq_model)[2])){
    to_remove = which(is.na(datSeq_model[,j])==TRUE)
    if(length(to_remove)>0){
      tmp_topPC=TopPC1[-to_remove,i]
      tmp_seq=datSeq_model[-to_remove,j]}
    else{
      tmp_topPC=TopPC1[,i]
      tmp_seq=datSeq_model[,j]  
    }
    mod=summary(lm(tmp_topPC~as.numeric(tmp_seq)))
    r2mat_seq[i,j]=mod$adj.r.squared
  }
}

pdf(file="Normalized_PC-heatmap.pdf")
pheatmap(r2mat_meta,clustering_method="average",cluster_rows=FALSE,display_numbers = TRUE,fontsize_number = 5)
pheatmap(t(r2mat_seq),clustering_method="average",fontsize = 6,cluster_cols=FALSE, display_numbers = TRUE,fontsize_number = 10)
dev.off()

########################## 
##### PC correction, for example, regress the first 5 PC, Batch, Sex etc. if applicable

Group <- as.numeric(as.factor(target$group))-1
#Sex=as.numeric(factor(target$Sex))
#Num.sort.x1k=as.numeric(factor(target$Num_Sort.x1k.))
#Batch=as.numeric(factor(target$Batch))
seqPC1 <- as.numeric(target$Seq.PC1)
seqPC2 <- as.numeric(target$Seq.PC2)
seqPC3 <- as.numeric(target$Seq.PC3)
seqPC4 <- as.numeric(target$Seq.PC4)
seqPC5 <- as.numeric(target$Seq.PC5)
seqPC6 <- as.numeric(target$Seq.PC6)

###PC1, PC4, and PC5 correlated with Seq PC, and PC2 and PC3 strongly correlated wth "Group"
modname="PC1-PC4-PC5"
X = model.matrix(~Group+seqPC1+seqPC4+seqPC5,data=target)
Y = datExpr

beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
b = as.data.frame(t(beta))
to_regress = (as.matrix(X[,3:dim(X)[2]]) %*% (as.matrix(beta[3:dim(X)[2],])))
datExprfit = datExpr - t(to_regress)

##########################
### MDS plot after PC correction

pdf(paste("MDSplot_PCAcorrection_",modname,".pdf",sep=""), width=10,height=14)
par(mfrow=c(2,2))
data <- datExprfit
ldat <- dist(t(data))
fit <- isoMDS(ldat, k=4)
x <- fit$points[,1]
y <- fit$points[,2]

for(i in 2:dim(conds)[2]){
  plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="VSD normalized genes with 2 counts for half of samples", type="p",pch=16,cex=1,col=clab[,i],xlim=c(min(x)*1.4,max(x)*1.4))
  text(x, y+(max(y)-min(y))/50, labels = row.names(t(data)), cex=sizeText, col=clab[,i])
  legend("topleft", legend=unique(as.factor(conds[,i])), fill=unique(clab[,i]),cex=sizeText)
}

dev.off()

##################
pdf(paste("MDSplot_PCAcorrection_",modname,"_Top1000.pdf",sep=""), width=10,height=14)
rowMadsT <- apply(datExprfit,1,mad)
ordT <- order(rowMadsT,decreasing=TRUE)
top500T <- head(ordT, n=500)

par(mfrow=c(2,2))
data <- datExpr[top500T,]
ldat <- dist(t(data))
fit <- isoMDS(ldat, k=4)
x <- fit$points[,1]
y <- fit$points[,2]

for(i in 2:dim(conds)[2]){
  plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="VSD normalized genes with 2 counts for half of samples", type="p",pch=16,cex=1,col=clab[,i],xlim=c(min(x)*1.4,max(x)*1.4))
  text(x, y+(max(y)-min(y))/50, labels = row.names(t(data)), cex=sizeText, col=clab[,i])
  legend("topleft", legend=unique(as.factor(conds[,i])), fill=unique(clab[,i]),cex=sizeText)
}

x <- fit$points[,3]
y <- fit$points[,4]

for(i in 2:dim(conds)[2]){
  plot(x, y, xlab="Coordinate 3", ylab="Coordinate 4", main="VSD normalized genes with 5 counts for half of samples", type="p",pch=16,cex=1,col=clab[,i],xlim=c(min(x)*1.4,max(x)*1.4))
  text(x, y+(max(y)-min(y))/50, labels = row.names(t(data)), cex=sizeText, col=clab[,i])
  legend("topleft", legend=unique(as.factor(conds[,i])), fill=unique(clab[,i]),cex=sizeText)
}

dev.off()

################################################

normExpr<-datExprfit
pdf(paste("PCA_logTMMNorm_PCA_corrected_",modname,".pdf",sep=""),height=20,width=24)

thisdat.HTSC <- t(scale(t(normExpr),scale=F))
PC.HTSC <- prcomp(thisdat.HTSC,center=F);
TopPC2 <- PC.HTSC$rotation[,1:5];
varexp <- (PC.HTSC$sdev)^2 / sum(PC.HTSC$sdev^2)
topvar <- varexp[1:5]
colnames(TopPC2) <- paste("Exp.PC\n",colnames(TopPC2)," (",signif(100*topvar[1:5],2),"%)",sep="")

pairs(cbind(TopPC2,pairsdat),col= cond,pch=19,upper.panel = panel.cor,main="Covariates and HTSeqCount Comparison -- |Spearman's rho| correlation values")
dev.off()

write.csv(datExprfit,paste("PCA_corrected_logNorm_exprn_",modname,".csv",sep=""))

r2mat_meta = matrix(NA,nrow=5,ncol=dim(conds)[2]-1)
r2mat_seq = matrix(NA,nrow=5,ncol=dim(target)[2]-dim(conds)[2])
rownames(r2mat_meta) <- rownames(r2mat_seq) <- colnames(TopPC2)
datMeta_model = target[,2:n]
datSeq_model=target[,-c(1:n)]

colnames(r2mat_meta) <- colnames(datMeta_model)
colnames(r2mat_seq) <- colnames(datSeq_model)

# Non-numeric factors vs seqQC
for(i in c(1:5)){
  for(j in c(1:dim(datMeta_model)[2])){
    to_remove = which(is.na(datMeta_model[,j])==TRUE)
    if(length(to_remove)>0){
      tmp_topPC=TopPC2[-to_remove,i]
      tmp_meta=datMeta_model[-to_remove,j]}
    else{
      tmp_topPC=TopPC2[,i]
      tmp_meta=datMeta_model[,j]  
    }
    mod_mat=model.matrix(~tmp_meta)[,-1]
    mod=summary(lm(tmp_topPC~mod_mat))
    r2mat_meta[i,j]=mod$adj.r.squared
  }
}

# numeric factors vs seqQC
for(i in c(1:5)){
  for(j in c(1:dim(datSeq_model)[2])){
    to_remove = which(is.na(datSeq_model[,j])==TRUE)
    if(length(to_remove)>0){
      tmp_topPC=TopPC2[-to_remove,i]
      tmp_seq=datSeq_model[-to_remove,j]}
    else{
      tmp_topPC=TopPC2[,i]
      tmp_seq=datSeq_model[,j]  
    }
    mod=summary(lm(tmp_topPC~as.numeric(tmp_seq)))
    r2mat_seq[i,j]=mod$adj.r.squared
  }
}

pdf(paste("Normalized_PC-heatmap_corrected_",modname,".pdf",sep=""))
pheatmap(r2mat_meta,clustering_method="average",cluster_rows=FALSE,display_numbers = TRUE,fontsize_number = 10)
pheatmap(t(r2mat_seq),clustering_method="average",fontsize = 6,cluster_cols=FALSE, display_numbers = TRUE,fontsize_number = 10)
dev.off()



