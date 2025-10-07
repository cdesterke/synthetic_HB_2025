list.files()

## load normalized RNAseq data
load("TNMVOOMfiltrated12386.rda")

## load transcriptome data
gse131<-read.csv("matrix.csv",h=T)
gse131<-gse131[complete.cases(gse131$gene),]


## gene symbol as row.names
library(transpipe15)
gse131<-filtermatrix(gse131)


## merge the two transcriptome datasets by row.names (gene symbols)
all<-merge(gse131,gse104,by="row.names")
colnames(all)
row.names(all)<-all$Row.names
all$Row.names<-NULL

## load phenotypes
pheno<-read.csv("pheno.csv",h=T,row.names=1)

boxplot(all)
all(row.names(pheno)==colnames(all))
data<-all[,row.names(pheno)]
all(row.names(pheno)==colnames(data))

## combat batch correction

library(sva)
batch = pheno$batch
mod = model.matrix(~1+as.factor(group), data=pheno)
edata = ComBat(dat=data,  batch=batch, mod=mod,par.prior=TRUE, prior.plots=TRUE)
boxplot(edata)


## quantiles normalization

library(preprocessCore)
norm_edata = normalize.quantiles(as.matrix(edata))
row.names(norm_edata)<-row.names(edata)
colnames(norm_edata)<-colnames(edata)

colramp = colorRampPalette(c(3,"darkblue",5))(8)
plot(density(norm_edata[,1]),col=colramp[1],lwd=3,ylim=c(0,.3))
for(i in 1:129){lines(density(norm_edata[,i]),lwd=3,col=colramp[i])}

boxplot(norm_edata)

save(norm_edata,file="quantile.rda")
write.csv(norm_edata,file="quantile.csv",row.names=T)


## transcriptome visualizations and DEG analysis (transpipe15)
pcatrans(edata,pheno,group="group",pal="Set1",alpha=1,names=F)

library(dplyr)
pheno%>%mutate(group2=ifelse(group=="normal","normal","tumor"))->pheno
res<-deg(norm_edata,pheno$group2,control="normal")
vollimma(res,nb=500,fc=1,p=0.01,size=4,alpha=1)

write.table(edata,file="superData.csv",row.names=T,sep="\t")

res%>%filter(adj.P.Val<=0.05 & abs(logFC)>=1)->sig

process<-reducedf(sig,norm_edata,n=1230)
pheno%>%select(group,batch)->annot
pcatrans(process,pheno,group="group",pal="Set1",alpha=1,names=F)
bestheat(process,annot,scale="row",font=8,rownames=F)


write.table(res,file="limmares.tsv",row.names=T,sep="\t")


## up regulated genes
sig%>%filter(logFC>0)->pos
process<-reducedf(pos,norm_edata,n=499)
pcatrans(process,pheno,group="group",pal="Set1",alpha=1,names=F)
bestheat(process,annot,scale="row",font=8,rownames=F)
write.table(pos,file="pos499.tsv",row.names=T,sep="\t")
#########################







