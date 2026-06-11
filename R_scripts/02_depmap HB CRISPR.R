# Données de dépendance génétique
#depmap_scores <- read_csv("CRISPR_gene_effect.csv")
load("depmap_scores.rda")
# Métadonnées des lignées
#sample_info <- read_csv("sample_info.csv")
load("depmap_samples.rda")

library(dplyr)
library(stringr)
# Extraire les hcc



lines_hcc <-sample_info %>% dplyr::filter(str_detect(Cellosaurus_NCIt_disease, regex("hepatocellular", ignore_case = TRUE)))
lines_hcc<-as.data.frame(lines_hcc)

id_hcc<-lines_hcc$depmap_id


depmap_hcc <- depmap_scores %>%
  dplyr::filter(depmap_id %in% id_hcc)



## hb
lines_hb <-sample_info %>% dplyr::filter(str_detect(Cellosaurus_NCIt_disease, regex("Hepatoblastoma", ignore_case = TRUE)))
lines_hb<-as.data.frame(lines_hb)


##id<-lines%>%  pull(depmap_id) %>% unique()
id_hb<-c("ACH-000739","ACH-000671","ACH-001021","ACH-002395")

# Filtrer les scores CRISPR pour ces lignées
# Garder uniquement les lignes correspondant aux lignées colorectales
depmap_hb <- depmap_scores %>%
  dplyr::filter(depmap_id %in% id_hb)



## cholangio
lines_ch <-sample_info %>% dplyr::filter(str_detect(Cellosaurus_NCIt_disease, regex("cholangiocarcinoma", ignore_case = TRUE)))
lines_ch<-as.data.frame(lines_ch)


##id<-lines%>%  pull(depmap_id) %>% unique()
id_ch<-lines_ch$depmap_id


# Filtrer les scores CRISPR pour ces lignées
# Garder uniquement les lignes correspondant aux lignées colorectales
depmap_ch <- depmap_scores %>%
  dplyr::filter(depmap_id %in% id_ch)









depmap_hcc$group<-"hepatocellular carcinoma"
depmap_hb$group<-"hepatoblastoma"
depmap_ch$group<-"cholangiocarcinoma"

load("New_SL_common73_genes_R1.rda")
depmap_all<-rbind(depmap_hb,depmap_hcc,depmap_ch)

depmap_all<-as.data.frame(depmap_all)
df<-read.csv("gene_centrality_42.csv",h=T)
vector<-df$Node
depmap_all<-depmap_all[depmap_all$gene_name%in%vector,]


ggplot(depmap_all,aes(cell_line,dependency))+geom_jitter(aes(color=group))+coord_flip()


depmap_all%>%filter(gene_name%in%common_genes)->depmap_sel

##tableau summary
library(dplyr)
library(tidyr)

stats_depmap <- depmap_sel %>%
  group_by(gene_name, group) %>%
  summarise(
    mean_dependency = mean(dependency, na.rm = TRUE),
    sd_dependency   = sd(dependency, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) 

#%>%  mutate(mean_sd = sprintf("%.3f ± %.3f", mean_dependency, sd_dependency))

stats_depmap
write.csv(stats_depmap,file="depmap73_summaryStats.csv",row.names=F)





library(dplyr)

stats <- depmap_all %>%
  group_by(cell_line) %>%
  summarise(
    mean_dep = mean(dependency, na.rm = TRUE),
    sd_dep   = sd(dependency, na.rm = TRUE)
  ) %>%
  arrange(mean_dep)   # tri du plus petit au plus grand

# Appliquer l’ordre aux données principales
depmap_all <- depmap_all %>%
  mutate(cell_line = factor(cell_line, levels = stats$cell_line))


library(ggplot2)

ggplot(depmap_all, aes(x = cell_line, y = dependency)) +
  
  # Jitter des valeurs individuelles
  geom_jitter(aes(color = group),
              width = 0.15, size = 2, alpha = 1) +
  
  # Barres mean ± SD
  geom_errorbar(
    data = stats,
    inherit.aes = FALSE,
    aes(x = cell_line,
        ymin = mean_dep - sd_dep,
        ymax = mean_dep + sd_dep),
    width = 0.25,
    color = "black",
    linewidth = 0.7
  ) +
  
  # Point de la moyenne
  geom_point(
    data = stats,
    inherit.aes = FALSE,
    aes(x = cell_line, y = mean_dep),
    color = "black",
    size = 3
  ) +
  
  coord_flip() +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 10),
    legend.position = "right"
  ) +
  labs(
    x = "Cell line (ranked by mean dependency)",
    y = "score CRISPR DepMap (Chronos/CERES)",
    title = "Dependency specificity comparison liver/biliary adult cancer"
  )




