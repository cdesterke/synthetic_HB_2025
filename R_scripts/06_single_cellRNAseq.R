library(Seurat)
load("harmony.rda")

x<-data[[]]
head(x)

input <- list(c("PLCG1","SORT1","PRKAA2","TSPAN5","ITGA6","PEG10","AXIN2","GJA5","EPCAM","PLK1","IGF2BP2","LTBP2","GPC3","ITGA2","LEF1","NOTUM"))
data<- AddModuleScore(
  object = data,
  features = input, name = 'score',
ctrl=16
)

library(ggplot2)
colgex=c("grey90",brewer.pal(7,"Reds"))
plotTheme=theme_classic(base_size=12)
FeaturePlot(data,reduction="umap",pt.size=0.1,features="en_score",min.cutoff = "q9",col=c("grey90","darkred"),split.by="sample_group")+
	plotTheme+coord_fixed()
Idents(data)<-"sample_group"
Cell.class

library(pals)
DimPlot(data,reduction="umap",group.by ="Cell.class_concise",pt.size=.1,cols=cols25(),label=F)

VlnPlot(data, features = c("en_score"), slot = "data", log = TRUE,pt.size=0.1,split.by="Cell.class_concise",col=cols25())
VlnPlot(data, features = c("en_score"), slot = "data", log = TRUE,pt.size=0.1,split.by="sample_group",col=cols25())

matrix<-data[["originalexp"]]$data

trans<-as.data.frame(t(matrix))

small<-trans[,colnames(trans)%in%genes_of_interest]

library(dplyr)


small%>%mutate(en_score=(PLCG1*0.429872123332352)+(SORT1*0.295592045676785)+(PRKAA2*0.195893492874045)
+(TSPAN5*0.195249942942498)+(ITGA6*0.195170798653957)+(PEG10*0.133784017006648)+
(AXIN2*0.132526109118132)+(GJA5*0.131870524080608)+(EPCAM*0.126018841613278)+
(PLK1*0.11118560539791)+(IGF2BP2*0.0452194527975264)+(LTBP2*0.0385934071695026)+
(GPC3*0.0337724056330341)+(ITGA2*0.0301561058236772)+(LEF1*0.0204097432488604)+
(NOTUM*0.00266070902875418))->small
small<-trans[,colnames(trans)%in%genes_of_interest]

data$en_score<-small$en_score
x<-data[[]]
# Liste des groupes uniques
groupes <- unique(x$sample_group)

# Créer un data frame pour stocker les résultats
resultats_df <- data.frame(
  Comparaison = character(),
  Moyenne_Groupe1 = numeric(),
  Moyenne_Groupe2 = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Boucle pour faire les tests t pour chaque paire de groupes
for (i in 1:(length(groupes) - 1)) {
  for (j in (i + 1):length(groupes)) {
    # Sélection des deux groupes à comparer
    groupe1 <- x$en_score[x$sample_group == groupes[i]]
    groupe2 <- x$en_score[x$sample_group == groupes[j]]
    
    # Test t
    test <- t.test(groupe1, groupe2)
    
    # Ajouter les résultats dans le data frame
    resultats_df <- rbind(resultats_df, data.frame(
      Comparaison = paste(groupes[i], "vs", groupes[j]),
      Moyenne_Groupe1 = mean(groupe1),
      Moyenne_Groupe2 = mean(groupe2),
      p_value = test$p.value
    ))
  }
}

# Afficher les résultats sous forme de data frame
print(resultats_df)
write.csv(resultats_df,file="comparisons_ENscore_between_groups.csv",row.names=F)





# Liste des groupes uniques cell types
groupes <- unique(x$Cell.class)

# Créer un data frame pour stocker les résultats
resultats_df <- data.frame(
  Comparaison = character(),
  Moyenne_Groupe1 = numeric(),
  Moyenne_Groupe2 = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Boucle pour faire les tests t pour chaque paire de groupes
for (i in 1:(length(groupes) - 1)) {
  for (j in (i + 1):length(groupes)) {
    # Sélection des deux groupes à comparer
    groupe1 <- x$en_score[x$Cell.class == groupes[i]]
    groupe2 <- x$en_score[x$Cell.class == groupes[j]]
    
    # Test t
    test <- t.test(groupe1, groupe2)
    
    # Ajouter les résultats dans le data frame
    resultats_df <- rbind(resultats_df, data.frame(
      Comparaison = paste(groupes[i], "vs", groupes[j]),
      Moyenne_Groupe1 = mean(groupe1),
      Moyenne_Groupe2 = mean(groupe2),
      p_value = test$p.value
    ))
  }
}

# Afficher les résultats sous forme de data frame
print(resultats_df)
write.csv(resultats_df,file="comparisons_ENscore_between_cell_classes.csv",row.names=F)


x$conca<-paste(x$sample_group,x$Cell.class,sep=".")


# Liste des groupes uniques cell types
groupes <- unique(x$conca)

# Créer un data frame pour stocker les résultats
resultats_df <- data.frame(
  Comparaison = character(),
  Moyenne_Groupe1 = numeric(),
  Moyenne_Groupe2 = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Boucle pour faire les tests t pour chaque paire de groupes
for (i in 1:(length(groupes) - 1)) {
  for (j in (i + 1):length(groupes)) {
    # Sélection des deux groupes à comparer
    groupe1 <- x$en_score[x$conca == groupes[i]]
    groupe2 <- x$en_score[x$conca == groupes[j]]
    
    # Test t
    test <- t.test(groupe1, groupe2)
    
    # Ajouter les résultats dans le data frame
    resultats_df <- rbind(resultats_df, data.frame(
      Comparaison = paste(groupes[i], "vs", groupes[j]),
      Moyenne_Groupe1 = mean(groupe1),
      Moyenne_Groupe2 = mean(groupe2),
      p_value = test$p.value
    ))
  }
}

# Afficher les résultats sous forme de data frame
print(resultats_df)
write.csv(resultats_df,file="comparisons_ENscore_between_conca.csv",row.names=F)





tumor<-resultats_df[grepl("Tumor",resultats_df$Comparaison),]

write.table(tumor,file="ElasticNet16tumor.tsv",row.names=F,sep="\t")




# Gènes d’intérêt
genes_of_interest <- c("PLCG1", "SORT1", "PRKAA2", "TSPAN5", "ITGA6", "PEG10",
                       "AXIN2", "GJA5", "EPCAM", "PLK1", "IGF2BP2", "LTBP2",
                       "GPC3", "ITGA2", "LEF1", "NOTUM")

# Vérifie quels gènes sont présents
genes_present <- genes_of_interest[genes_of_interest %in% rownames(data)]

# Extraction des données normalisées (log-normalisées par défaut dans le slot "data")

normalized_expr <- GetAssayData(data, assay = "originalexp", slot = "data")[genes_present, ]
# Aperçu
normalized_expr[1:5, 1:5]
##

library(ggplot2)

library(RColorBrewer)

colgex=c("grey90",brewer.pal(7,"Reds"))
plotTheme=theme_classic(base_size=16)


DotPlot(data,group.by="sample_group",features=genes_of_interest)+
coord_flip() + scale_color_gradientn(colors = colgex) +
  theme(axis.text.x = element_text(angle = -90, vjust = 0))

library(dplyr)
x%>%select(1:10)->x

df<-as.data.frame(t(normalized_expr))

all<-merge(x,small,by="row.names")

all(row.names(x)==row.names(small))

all$conca<-x$conca

write.csv(all,file="dataDL.csv",row.names=F)

save(all,file="dataDL.rda")



## select tumor and hepa from conca

all%>%filter(grepl("background.Hepatocyte|tumor.Tumor cell",conca))->df

df%>%select(conca,all_of(genes_of_interest))->df2
table(df2$conca)


library(dplyr)

library(dplyr)

df2 <- df2 %>%
  mutate(conca_bin = recode(conca,
                            "tumor.Tumor cell" = 1,
                            "background.Hepatocyte" = 0,
                            .default = NA_real_))


df2$conca<-NULL

df2%>%relocate(conca_bin)->df


write.csv(df,file="dataset.csv",row.names=F)





