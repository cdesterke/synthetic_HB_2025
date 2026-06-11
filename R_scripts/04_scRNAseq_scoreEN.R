library(Seurat)
load("harmony.rda")

x<-data[[]]
head(x)

genes_of_interest<-c("PEG10","CRIM1","SMYD2","FABP4","FSD1L","PEG3","PLCB4","RHOBTB1","UTRN","KIT","HDAC11","SEMA7A","GREB1","ANKRD50","ZNF233","CORO2A")




input <- list(c("PEG10","CRIM1","SMYD2","FABP4","FSD1L","PEG3","PLCB4","RHOBTB1","UTRN","KIT","HDAC11","SEMA7A","GREB1","ANKRD50","ZNF233","CORO2A"))
data<- AddModuleScore(
  object = data,
  features = input, name = 'en_score',
ctrl=16
)

library(RColorBrewer)
library(ggplot2)
colgex=c("grey90",brewer.pal(7,"Reds"))
plotTheme=theme_classic(base_size=12)
FeaturePlot(data,reduction="umap",pt.size=0.1,features="en_score1",min.cutoff = "q9",col=c("grey90","darkred"),split.by="sample_group")+
	plotTheme+coord_fixed()
Idents(data)<-"sample_group"
Cell.class

library(pals)
DimPlot(data,reduction="umap",group.by ="Cell.class_concise",pt.size=.1,cols=cols25(),label=F)

VlnPlot(data, features = c("en_score1"), slot = "data", log = TRUE,pt.size=0.1,split.by="Cell.class_concise",col=cols25())
VlnPlot(data, features = c("en_score1"), slot = "data", log = TRUE,pt.size=0.1,split.by="sample_group",col=cols25())

matrix<-data[["originalexp"]]$data

trans<-as.data.frame(t(matrix))

small<-trans[,colnames(trans)%in%genes_of_interest]

library(dplyr)





small<-trans[,colnames(trans)%in%genes_of_interest]

data$en_score<-small$en_score
x<-data[[]]
# Liste des groupes uniques
groupes <- unique(all$sample_group)

# Créer un data frame pour stocker les résultats
resultats_df <- data.frame(
  Comparaison = character(),
  Median_Groupe1 = numeric(),
  Median_Groupe2 = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Boucle pour faire les tests t pour chaque paire de groupes
resultats_df <- data.frame()

for (i in 1:(length(groupes) - 1)) {
  for (j in (i + 1):length(groupes)) {
    
    # Sélection des deux groupes à comparer
    groupe1 <- all$en_score[all$sample_group == groupes[i]]
    groupe2 <- all$en_score[all$sample_group == groupes[j]]
    
    # Test Wilcoxon (non paramétrique)
    test <- wilcox.test(groupe1, groupe2)
    
    # Ajouter les résultats dans le data frame
    resultats_df <- rbind(resultats_df, data.frame(
      Comparaison = paste(groupes[i], "vs", groupes[j]),
      Median_Groupe1 = median(groupe1, na.rm = TRUE),
      Median_Groupe2 = median(groupe2, na.rm = TRUE),
      p_value = test$p.value
    ))
  }
}


# Afficher les résultats sous forme de data frame
print(resultats_df)
write.csv(resultats_df,file="NEW_comparisons_ENscore_between_groups.csv",row.names=F)





# Liste des groupes uniques cell types
groupes <- unique(all$Cell.class)

# Créer un data frame pour stocker les résultats
resultats_df <- data.frame(
  Comparaison = character(),
  Median_Groupe1 = numeric(),
  Median_Groupe2 = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Boucle pour faire les tests t pour chaque paire de groupes
for (i in 1:(length(groupes) - 1)) {
  for (j in (i + 1):length(groupes)) {
    # Sélection des deux groupes à comparer
    groupe1 <- all$en_score[all$Cell.class == groupes[i]]
    groupe2 <- all$en_score[all$Cell.class == groupes[j]]
    
    # Test t
    test <- wilcox.test(groupe1, groupe2)
    
    # Ajouter les résultats dans le data frame
    resultats_df <- rbind(resultats_df, data.frame(
      Comparaison = paste(groupes[i], "vs", groupes[j]),
      Median_Groupe1 = median(groupe1),
      Median_Groupe2 = median(groupe2),
      p_value = test$p.value
    ))
  }
}

# Afficher les résultats sous forme de data frame
print(resultats_df)
write.csv(resultats_df,file="NEW_comparisons_ENscore_between_cell_classes.csv",row.names=F)


all$conca<-paste(all$sample_group,all$Cell.class,sep=".")


# Liste des groupes uniques cell types
groupes <- unique(all$conca)

# Créer un data frame pour stocker les résultats
resultats_df <- data.frame(
  Comparaison = character(),
  Median_Groupe1 = numeric(),
  Median_Groupe2 = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Boucle pour faire les tests t pour chaque paire de groupes
for (i in 1:(length(groupes) - 1)) {
  for (j in (i + 1):length(groupes)) {
    # Sélection des deux groupes à comparer
    groupe1 <- all$en_score[all$conca == groupes[i]]
    groupe2 <- all$en_score[all$conca == groupes[j]]
    
    # Test t
    test <- wilcox.test(groupe1, groupe2)
    
    # Ajouter les résultats dans le data frame
    resultats_df <- rbind(resultats_df, data.frame(
      Comparaison = paste(groupes[i], "vs", groupes[j]),
      Median_Groupe1 = mean(groupe1),
      Median_Groupe2 = mean(groupe2),
      p_value = test$p.value
    ))
  }
}

# Afficher les résultats sous forme de data frame
print(resultats_df)
write.csv(resultats_df,file="NEWcomparisons_ENscore_between_conca.csv",row.names=F)





tumor<-resultats_df[grepl("Tumor",resultats_df$Comparaison),]

write.table(tumor,file="ElasticNet16tumor.tsv",row.names=F,sep="\t")




# Gènes d’intérêt


# Vérifie quels gènes sont présents
genes_present <- genes_of_interest[genes_of_interest %in% rownames(data)]

# Extraction des données normalisées (log-normalisées par défaut dans le slot "data")

normalized_expr <- GetAssayData(data, assay = "originalexp", layer = "data")[genes_present, ]
# Aperçu
normalized_expr[1:5, 1:5]
##

library(ggplot2)

library(RColorBrewer)

colgex=c("grey90",brewer.pal(7,"Reds"))
plotTheme=theme_classic(base_size=16)


DotPlot(data,group.by="Cell.class_concise",features=genes_of_interest)+
coord_flip() + scale_color_gradientn(colors = colgex) +
  theme(axis.text.x = element_text(angle = -90, vjust = 0))
x<-data[[]]
library(dplyr)
x%>%select(1:10)->small
small$en_score<-x$en_score1
df<-as.data.frame(t(normalized_expr))
all(row.names(df)==row.names(small))
all<-merge(df,small,by="row.names")



all$conca<-x$conca

write.csv(all,file="data16genesNew.csv",row.names=F)

save(all,file="data16genesNew.rda")


















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


write.csv(df,file="datasetNEW.csv",row.names=F)



