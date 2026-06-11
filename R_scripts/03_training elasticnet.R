list.files()
load("New_SL_common73_genes_R1.rda")
load("ESET_GSE104766.rda")
small<-data[row.names(data)%in%common_genes,]
all(colnames(small)==row.names(pheno))

trans<-t(small)
scaled<-scale(trans)

apply(scaled, 2, mean) |> head()
apply(scaled, 2, sd) |> head()


library(glmnet)
y <- as.factor(pheno$tissue)
y<-relevel(y,ref="Normal liver")
X<-as.matrix(scaled)

library(caret)
trainIndex <- createDataPartition(y, p = 0.65, list = FALSE)
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
library(glmnet)
library(dplyr)
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
library(pals)

results <- results %>%
  mutate(log_lambda = log10(lambda))

ggplot(results, aes(x = log_lambda, y = auc, color = as.factor(alpha))) +
  geom_line() +
  geom_point() +
  
  # ligne verticale sur le meilleur lambda
  geom_vline(
    xintercept = log10(best_lambda),
    linetype = "dashed",
    color = "red",
    size = 1
  ) +
  
  # annotation texte
  annotate(
    "text",
    x = log10(best_lambda),
    y = max(results$auc),
    label = paste0(
      "Best alpha = ", best_alpha,
      "\nBest lambda = ", signif(best_lambda, 3),
      "\nBest AUC = ", round(best_auc, 3)
    ),
    hjust = -1,
    vjust = 1,
    size = 5,
    color = "red",
    fontface = "bold"
  ) +
  
  scale_color_manual(values = glasbey()) +
  labs(
    title = "GSE104766 internal split 0.65/0.35",
    x = "log10(lambda)",
    y = "AUC",
    color = "Alpha"
  ) +
  theme_minimal() +
  theme(text = element_text(size = 16))



### elastic net 0.7 alpha 

fit <- glmnet(X_train, y_train, family = "binomial",alpha=0.7)
plot(fit)

cvfit <- cv.glmnet(X_train, y_train, family = "binomial",alpha=0.7)
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

write.table(selcoef,file="selcoef16New.tsv",row.names=T,sep="\t")



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




## test 16 genes

all<-cbind(y_test,X_test)
all<-as.data.frame(all)
all%>%mutate(score=(PEG10*0.443076288090222)+(CRIM1*0.431377193400511)+(SMYD2*0.420443438316175)+
(FABP4*0.300037633307907)+(FSD1L*0.295893162984562)+(PEG3*0.290662361320897)+
(PLCB4*0.286313859384664)+(RHOBTB1*0.235763430544155)+(UTRN*0.136909964138547)+
(KIT*0.103070748953849)+(HDAC11*0.0925106341892726)+(SEMA7A*0.0869042849504362)+
(GREB1*0.0770673190189149)+(ANKRD50*0.0422888470086895)+(ZNF233*0.0237716008330501)+
(CORO2A*0.000218362627059983))->all

library(cutpointr)


cp <- cutpointr(all, score,  y_test,method = maximize_metric, metric = sum_sens_spec)

plot(cp)
cp
all%>%mutate(en_cat=ifelse(score>=0.153905,"HIGH","low"))->all

chisq.test(table(all$en_cat,all$y_test))


## mosaicplot

library(vcd)
all$y_test<-as.character(all$y_test)
all <- all %>%
  mutate(
    y_valid = case_when(
      y_test == "2" ~ "HB",
      y_test == "1" ~ "NL",
      TRUE ~ NA_character_
    )
  )

struct <- structable(~ en_cat+y_valid,data = all)
mosaic(struct, , direction = "h", pop = FALSE,colorize = T, shade = TRUE)
       #gp = gpar(fill = matrix(c("red","grey90" , "grey90","grey90" , "grey90", "green3"), 2, 3)))
labeling_cells(text = as.table(struct), margin = 0)(as.table(struct))
chisq.test(struct)


## Boxplot 
library(ggplot2)
library(ggbeeswarm)
all%>%dplyr::rename(en_score="score")->all

all<-all[,1:77]

ggplot(all,aes(y_valid,en_score))+geom_boxplot(outlier.shape=NA) + 
  scale_fill_brewer(palette="Set1")+
  geom_point(aes(fill=factor(y_valid),size=1),shape = 21, alpha = 1, position = position_dodge2(width = .5))+
  theme_classic(base_size=16) +
  theme(legend.position = "none")+xlab("liver tissues")+ylab("Elastic-net score")+ggtitle("GSE104766 validation split 0.35")




library(ggplot2)
library(ggpubr)

ggplot(all, aes(x = y_valid, y = en_score, fill = y_valid)) +
  
  # Boxplot propre
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.9, color = "black") +
  
  # Points jitter
  geom_jitter(
    shape = 21, size = 4, alpha = 1,
    width = 0.15, stroke = 0.3
  ) +
  
  # Palette
  scale_fill_brewer(palette = "Set2") +
  
  # Test t + barre + p-value formatée
  stat_compare_means(
    method = "t.test",
    label = "p.format",
    size = 5,
    bracket.size = 0.9,
    tip.length = 0.02,
    label.y = max(all$en_score, na.rm = TRUE) * 1.05
  ) +
  
  labs(
    x = "Liver tissues",
    y = "Elastic-net score",
    title = "GSE104766 – validation split 0.35"
  ) +
  
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )


