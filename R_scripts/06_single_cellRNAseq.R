library(Seurat)
load("harmony.rda")

x<-data[[]]
head(x)

input <- list(c("PLCG1","SORT1","PRKAA2","TSPAN5","ITGA6","PEG10","AXIN2","GJA5","EPCAM","PLK1","IGF2BP2","LTBP2","GPC3","ITGA2","LEF1","NOTUM"))
data<- AddModuleScore(
  object = data,
  features = input, name = 'en_score',
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


x<-data[[]]
# list of groups to compare
groupes <- unique(x$sample_group)

# create dataframe
resultats_df <- data.frame(
  Comparaison = character(),
  Moyenne_Groupe1 = numeric(),
  Moyenne_Groupe2 = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# loop of ttest
for (i in 1:(length(groupes) - 1)) {
  for (j in (i + 1):length(groupes)) {
    # Sélection des deux groupes à comparer
    groupe1 <- x$en_score[x$sample_group == groupes[i]]
    groupe2 <- x$en_score[x$sample_group == groupes[j]]
    
    # Test t
    test <- t.test(groupe1, groupe2)
    
    # add data to df
    resultats_df <- rbind(resultats_df, data.frame(
      Comparaison = paste(groupes[i], "vs", groupes[j]),
      Moyenne_Groupe1 = mean(groupe1),
      Moyenne_Groupe2 = mean(groupe2),
      p_value = test$p.value
    ))
  }
}


write.csv(resultats_df,file="comparisons_ENscore_between_groups.csv",row.names=F)





# list of cell types
groupes <- unique(x$Cell.class)

# creaste df
resultats_df <- data.frame(
  Comparaison = character(),
  Moyenne_Groupe1 = numeric(),
  Moyenne_Groupe2 = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# loop of ttest between cell types
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


write.csv(resultats_df,file="comparisons_ENscore_between_cell_classes.csv",row.names=F)


x$conca<-paste(x$sample_group,x$Cell.class,sep=".")


# list of groups
groupes <- unique(x$conca)

# create df
resultats_df <- data.frame(
  Comparaison = character(),
  Moyenne_Groupe1 = numeric(),
  Moyenne_Groupe2 = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# loop of ttests
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


write.csv(resultats_df,file="comparisons_ENscore_between_conca.csv",row.names=F)


## filter only tumor test
tumor<-resultats_df[grepl("Tumor",resultats_df$Comparaison),]
write.table(tumor,file="ElasticNet16tumor.tsv",row.names=F,sep="\t")




# genes of interest
genes_of_interest <- c("PLCG1", "SORT1", "PRKAA2", "TSPAN5", "ITGA6", "PEG10",
                       "AXIN2", "GJA5", "EPCAM", "PLK1", "IGF2BP2", "LTBP2",
                       "GPC3", "ITGA2", "LEF1", "NOTUM")

genes_present <- genes_of_interest[genes_of_interest %in% rownames(data)]

# extract single cell expression
normalized_expr <- GetAssayData(data, assay = "originalexp", slot = "data")[genes_present, ]

library(ggplot2)
library(RColorBrewer)

colgex=c("grey90",brewer.pal(7,"Reds"))
plotTheme=theme_classic(base_size=16)
DotPlot(data,group.by="sample_group",features=genes_of_interest)+
coord_flip() + scale_color_gradientn(colors = colgex) +
  theme(axis.text.x = element_text(angle = -90, vjust = 0))

library(dplyr)
x%>%select(1:10)->x

## extract data for deeplearning
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
df2 <- df2 %>%
  mutate(conca_bin = recode(conca,
                            "tumor.Tumor cell" = 1,
                            "background.Hepatocyte" = 0,
                            .default = NA_real_))

df2$conca<-NULL
df2%>%relocate(conca_bin)->df

write.csv(df,file="dataset.csv",row.names=F)






