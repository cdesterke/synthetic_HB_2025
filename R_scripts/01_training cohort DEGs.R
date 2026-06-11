list.files()

load("quantile.rda")
colramp = colorRampPalette(c(3,"darkblue",5))(24)
plot(density(norm[,1]),col=colramp[1],lwd=3,ylim=c(0,.2))
	for(i in 1:ncol(norm)){lines(density(norm[,i]),lwd=1,col=colramp[i])}

pheno<-read.csv("pheno.csv",h=T,row.names=1)

data<-norm

all(colnames(data)==row.names(pheno))

data<-data[,row.names(pheno)]
library(dplyr)
## remove cell lines
pheno%>%filter(tissue!="Hepatoblastoma cell line")->pheno
data<-data[,row.names(pheno)]

save(data,pheno,file="ESET_GSE104766.rda")


library(transpipe15)
load("vector.rda")
dim(data)
datasel<-data[row.names(data)%in%vector,]

res<-deg(datasel,pheno$tissue,control="Normal liver")



## filter significant DEGs
sig<-filtresig(res)

## draw volcanoplot
vollimma(res,nb=500,fc=1,p=0.01,size=4,alpha=1)

sig%>%filter(logFC>=1)->pos

## subset matrix to significant genes
process<-reducedf(pos,datasel,n=789)


pca_plot(process,
         pheno,
         group = "tissue",
         point_shape = "tissue",
         pal = "Set2",
         ellipse_conf = 0.75,
         point_size = 5,
         base_size = 18,alpha=1,
         title = "GSE104766/789 genes up regulated in tumors",
         show_manova = TRUE)

pheno%>%select(tissue)->annot
bestheat(process,annot,rownames=F,font=12)


save(pos,file="pos789_training.rda")
