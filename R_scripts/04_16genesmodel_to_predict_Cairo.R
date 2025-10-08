
pheno<-read.table("pheno_tumor.tsv",h=T,sep="\t")



data<-read.table("resultsScores.csv",h=T,sep=",",row.names=1)





all<-merge(data,pheno,bdatay="row.names")


library(dplyr)
all%>%select(Row.names,en_score,en_cat,tissue.x,group.x,dataset,chic_risk_stratification,
		clinical_course,ctnnbi_gene_alteration.ch1,metastasis,metastasis_code,age_months,event,
		event_code,histological_type,gender,pretext_stage,cairo,cairo_code)->df


## multivariate model
model<-glm(cairo_code~en_score+histological_type+age_months+gender+pretext_stage,data=df,family="binomial")
summary(model)

library(broom)
library(broom.helpers)
library(GGally)

ggcoef_model(model,exponentiate=T,colour_guide=TRUE)+
scale_color_brewer(palette="Set2")+theme(text=element_text(size=18))+theme(legend.position="none")


# Coefficient for the predictor
coef_predictor <- coef(model)["en_score"]

# Odds Ratio
odds_ratio <- exp(coef_predictor)
odds_ratio

library(regplot)
regplot(model, clickable=F, points=T, droplines=T,rank="sd") 



genes_of_interest <- c("PLCG1", "SORT1", "PRKAA2.y", "TSPAN5", "ITGA6", "PEG10",
                       "AXIN2", "GJA5", "EPCAM", "PLK1", "IGF2BP2", "LTBP2",
                       "GPC3", "ITGA2", "LEF1", "NOTUM")


inter<-intersect(genes_of_interest,colnames(all))

all%>%select(Row.names,all_of(inter))->small
row.names(small)<-small$Row.names
small$Row.names<-NULL
small%>%dplyr::rename(PRKAA2="PRKAA2.y")->small
trans<-as.data.frame(t(small))
pretext_stage,metastasis
all%>%select(cairo,gender,)->pheno

row.names(pheno)<-colnames(trans)

library(transpipe15)

bestheat(trans,pheno,font=10,scale="row")

pcatransellipses(trans,pheno,group="cairo",alpha=1,x=1,y=2,level=0.5)


df<-cbind(small,pheno)

df$cairo<-as.factor(df$cairo)