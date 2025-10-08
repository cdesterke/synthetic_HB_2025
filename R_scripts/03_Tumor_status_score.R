list.files()

load("quantile.rda")

pheno<-read.csv("pheno.csv",h=T,row.names=1)
library(dplyr)
pheno%>%filter(group!="cell.line")->pheno

data<-norm_edata[,row.names(pheno)]


sig<-read.csv("gene_centrality_42.csv",h=T)

input<-sig$Node

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
y <- as.factor(pheno$group)
X<-as.matrix(trans)

library(caret)
trainIndex <- createDataPartition(y, p = 0.6, list = FALSE)
X_train <- X[trainIndex,]
y_train <- y[trainIndex]
X_test <- X[-trainIndex,]
y_test <- y[-trainIndex]


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



### elastic net 0.7 alpha 

fit <- glmnet(X, y, family = "binomial",alpha=0.7)
plot(fit)

cvfit <- cv.glmnet(X, y, family = "binomial",alpha=0.7)
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

write.table(selcoef,file="selcoef.tsv",row.names=T,sep="\t")



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


(PLCG1*0.429872123332352)+(SORT1*0.295592045676785)+(PRKAA2*0.195893492874045)
+(TSPAN5*0.195249942942498)+(ITGA6*0.195170798653957)+(PEG10*0.133784017006648)+
(AXIN2*0.132526109118132)+(GJA5*0.131870524080608)+(EPCAM*0.126018841613278)+
(PLK1*0.11118560539791)+(IGF2BP2*0.0452194527975264)+(LTBP2*0.0385934071695026)+
(GPC3*0.0337724056330341)+(ITGA2*0.0301561058236772)+(LEF1*0.0204097432488604)+
(NOTUM*0.00266070902875418)
 
 


 


all%>%mutate(en_score=(PLCG1*0.429872123332352)+(SORT1*0.295592045676785)+(PRKAA2*0.195893492874045)
+(TSPAN5*0.195249942942498)+(ITGA6*0.195170798653957)+(PEG10*0.133784017006648)+
(AXIN2*0.132526109118132)+(GJA5*0.131870524080608)+(EPCAM*0.126018841613278)+
(PLK1*0.11118560539791)+(IGF2BP2*0.0452194527975264)+(LTBP2*0.0385934071695026)+
(GPC3*0.0337724056330341)+(ITGA2*0.0301561058236772)+(LEF1*0.0204097432488604)+
(NOTUM*0.00266070902875418))->all




library(cutpointr)


cp <- cutpointr(all, en_score,  group,method = maximize_metric, metric = sum_sens_spec)

plot(cp)
cp
all%>%mutate(en_cat=ifelse(en_score>=13.2136,"HIGH","low"))->all

chisq.test(table(all$en_cat,all$group))

## mosaicplot

library(vcd)

struct <- structable(~ en_cat+group,data = all)
mosaic(struct, , direction = "h", pop = FALSE,colorize = T, shade = TRUE)
       #gp = gpar(fill = matrix(c("red","grey90" , "grey90","grey90" , "grey90", "green3"), 2, 3)))
labeling_cells(text = as.table(struct), margin = 0)(as.table(struct))
chisq.test(struct)



## Boxplot 
library(ggplot2)
library(ggbeeswarm)

ggplot(all,aes(group,en_score))+geom_boxplot(outlier.shape=NA) + 
  scale_fill_brewer(palette="Set1")+
  geom_point(aes(fill=factor(group),size=1),shape = 21, alpha = 1, position = position_dodge2(width = .5))+
  theme_classic(base_size=16) +
  theme(legend.position = "none")+xlab("liver tissues")+ylab("Elastic-net score")+ggtitle("hepatoblastoma expression")





final<-merge(all,pheno,by="row.names")

str(final)


write.csv(final,file="resultsScores.csv",row.names=F)


library(multirocauc)

all$group<-as.factor(all$group)

geneset<-row.names(pos)


save(df,file="df.rda")
df<-trans[,colnames(trans)%in%geneset]
all(row.names(df)==row.names(pheno))
df$group<-pheno$group
df$group<-as.factor(df$group)


roc.list<-roclist(df,geneset,outcome="group")
roc.list

rocplot(roc.list,line=1,title="Hepatoblastoma tumor status",police=12)


save(final,file="final.rda")
