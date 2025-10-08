
pheno<-read.table("pheno_tumor.tsv",h=T,sep="\t")

pheno%>%select(42:58)->annot

data<-read.table("resultsScores.csv",h=T,sep=",",row.names=1)





all<-merge(data,pheno,by="row.names")

small<-all[,2:43]
row.names(small)<-all$Row.names
X<-as.matrix(small)

df<-data[row.names(data)%in%input,]

trans<-as.data.frame(t(df))

all<-cbind(trans,pheno$group)
geneset<-colnames(trans)
library(multirocauc)
all%>%dplyr::rename(group="pheno$group")->all

all$group<-as.factor(all$group)
roc.list<-roclist(all,geneset,outcome="group")
roc.list

rocplot(roc.list,line=1,title="Hepatoblastoma tumor status",police=8)

library(glmnet)
y <- as.factor(all$metastasis)
X<-as.matrix(trans)

library(caret)
trainIndex <- createDataPartition(y, p = 0.65, list = FALSE)
X_train <- X[trainIndex,]
y_train <- y[trainIndex]
X_test <- X[-trainIndex,]
y_test <- y[-trainIndex]

library(glmnet)
# sequence for tuning of lambda and alpha
lambda_grid <- 10^seq(3, -2, by = -0.1)
alpha_grid <- seq(0.1, 0.9, by = 0.1)

# initialisation of variables
results <- expand.grid(alpha = alpha_grid, lambda = lambda_grid)
results$auc <- NA
library(pROC)
# loop for lambda and alpha tuning
for (i in 1:nrow(results)) {
  alpha <- results$alpha[i]
  lambda <- results$lambda[i]
  
  model <- glmnet(X_train, y_train, alpha = alpha, lambda = lambda, family = "binomial")
  
  # Prédictions on test data
  probs <- predict(model, s = lambda, newx = X_test, type = "response")
  
  # compute AUC
  roc_obj <- roc(y_test, as.vector(probs))
  auc <- auc(roc_obj)
  
  # data results
  results$auc[i] <- auc
}


# Selection of best parameters
best_params <- results[which.max(results$auc),]
best_params
best_alpha <- best_params$alpha

best_lambda <- best_params$lambda
best_auc <- best_params$auc

library(ggplot2)
# tuning visualisation
results <- results %>%
  mutate(log_lambda = log10(lambda))
library(pals)
ggplot(results, aes(x = log_lambda, y = auc, color = as.factor(alpha))) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = log10(best_lambda), linetype = "dashed", color = "red",size=1) +
	scale_color_manual(values=glasbey())+
  labs(title = "AUC depending on lambda and alpha",
       x = "log10(lambda)",
       y = "AUC",
       color = "Alpha") +
  theme_minimal()+theme(text = element_text(size = 16))



### elastic net 0.1 alpha 

fit <- glmnet(X, y, family = "binomial",alpha=0.1)
plot(fit)

cvfit <- cv.glmnet(X, y, family = "binomial",alpha=0.1)
plot(cvfit)
coefficients<-coef(fit, s = cvfit$lambda.min)
coefficients<-as.matrix(coefficients)
coefficients<-as.data.frame(coefficients)

library(dplyr)
colnames(coefficients)<-"coef"
coefficients ->selcoef

selcoef$gene<-row.names(selcoef)
selcoef%>%arrange(desc(coef))->selcoef
selcoef%>%filter(gene != "(Intercept)")->selcoef
selcoef
selcoef%>%filter(coef>0)->pos
pos
dim(pos)

write.table(selcoef,file="selcoefmetastasis.tsv",row.names=T,sep="\t")



library(ggplot2)

ggplot(data = pos, aes(x = reorder(gene, coef), y = coef, fill = gene)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  xlab("Genes") +
  ylab("Elasticnet coefficients") +
  geom_text(aes(label = round(coef, 3)), hjust = 0, vjust = 0.5, color = "darkblue",
            position = position_dodge(0), size = 4, angle = 0) +
  scale_fill_manual(values = rainbow(length(unique(pos$gene)))) +
  ggtitle("") +
  theme(text = element_text(size = 14),
        legend.position = "none")




id<-pos$gene
beta<-pos$coef

equation<-paste(id,beta,sep="*",collapse=")+(")
equation<-paste("(",equation,")",sep="")
equation


(MDK*0.086526467052398)+(EPS8L3*0.0473156316066736)+(PLK1*0.0387226254171079)+
(PLCG1*0.0255997665905855)+(GSTP1*0.0235128091238549)+(PYCR1*0.0202066021506857)+
(FOXM1*0.0184254570167794)+(SLC2A1*0.0113882786355533)+(APLN*0.00165711286195314)+
(DLK1*0.000793685290516393)


 
 


 


all%>%mutate(meta10_score=(MDK*0.086526467052398)+(EPS8L3*0.0473156316066736)+(PLK1*0.0387226254171079)+
(PLCG1*0.0255997665905855)+(GSTP1*0.0235128091238549)+(PYCR1*0.0202066021506857)+
(FOXM1*0.0184254570167794)+(SLC2A1*0.0113882786355533)+(APLN*0.00165711286195314)+
(DLK1*0.000793685290516393))->all




library(cutpointr)


cp <- cutpointr(all, meta10_score,  metastasis,method = maximize_metric, metric = sum_sens_spec)

plot(cp)
cp
all%>%mutate(meta_cat=ifelse(meta10_score>=1.93086,"HIGH","low"))->all

chisq.test(table(all$meta_cat,all$metastasis))

## mosaicplot

library(vcd)

struct <- structable(~ meta_cat+metastasis,data = all)
mosaic(struct, , direction = "h", pop = FALSE,colorize = T, shade = TRUE)
       #gp = gpar(fill = matrix(c("red","grey90" , "grey90","grey90" , "grey90", "green3"), 2, 3)))
labeling_cells(text = as.table(struct), margin = 0)(as.table(struct))
chisq.test(struct)



## Boxplot 
library(ggplot2)
library(ggbeeswarm)

ggplot(all,aes(metastasis,meta10_score))+geom_boxplot(outlier.shape=NA) + 
  scale_fill_brewer(palette="Set1")+
  geom_point(aes(fill=factor(gender),size=1),shape = 21, alpha = 1, position = position_dodge2(width = .5))+
  theme_classic(base_size=16) +
  theme(legend.position = "right")+xlab("metastasis status")+ylab("metastasis score")+ggtitle("")








write.csv(all,file="metastisis_Scores.csv",row.names=F)


library(multirocauc)

all$metastasis<-as.factor(all$metastasis)

geneset<-row.names(pos)



df<-all[,colnames(all)%in%geneset]
all(row.names(df)==row.names(pheno))
df$metastasis<-all$metastasis
df$metastasis<-as.numeric(df$metastasis)

## multivariate model
model<-glm(metastasis~meta10_score+histological_type+age_months+gender+pretext_stage,data=all,family="binomial")
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

library(broom)
export <- tidy(model, exponentiate = TRUE, conf.int = TRUE)
write.table(export,file="modelCI95metastasis.tsv",sep="\t",row.names=F)

genes_of_interest <- row.names(pos)


inter<-intersect(genes_of_interest,colnames(all))

all%>%select(Row.names,all_of(inter))->small
row.names(small)<-small$Row.names
small$Row.names<-NULL
small%>%dplyr::rename(PRKAA2="PRKAA2.y")->small
trans<-as.data.frame(t(small))
pretext_stage,metastasis
all%>%select(metastasis,gender)->pheno

row.names(pheno)<-colnames(trans)

library(transpipe15)

bestheat(trans,pheno,font=10,scale="row")

pcatrans(trans,pheno,group="metastasis",alpha=1)


df<-cbind(small,pheno)

df$cairo<-as.factor(df$cairo)



save(final,file="final.rda")


library(Publish)
all<-read.csv("metastisis_Scores.csv",h=T)
all$meta_cat<-as.factor(all$meta_cat)
all$meta_cat<-relevel(all$meta_cat,ref="low")
u<-univariateTable(meta_cat~age_months+chic_risk_stratification+histological_type+gender+pretext_stage+cairo+metastasis,data=all)
results<-summary(u)
write.table(results,file="PublishMeta.tsv",row.names=F,sep="\t")

save(model,file="modelMeta.rda")








