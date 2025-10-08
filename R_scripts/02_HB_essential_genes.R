# Données de dépendance génétique
#depmap_scores <- read_csv("CRISPR_gene_effect.csv")
load("depmap_scores.rda")
# Métadonnées des lignées
#sample_info <- read_csv("sample_info.csv")
load("depmap_samples.rda")



# Extraire les hb

library(stringr)
lines <-sample_info %>% filter(str_detect(Cellosaurus_NCIt_disease, regex("Hepatoblastoma", ignore_case = TRUE)))
lines<-as.data.frame(lines)
write.csv(sample_info,file="metadata.csv",row.names=F)

##id<-lines%>%  pull(depmap_id) %>% unique()
id<-c("ACH-000739","ACH-000671")

# Filtrer les scores CRISPR pour ces lignées
# Garder uniquement les lignes correspondant aux lignées colorectales
depmap <- depmap_scores %>%
  filter(depmap_id %in% id)



## identifier les genes essentiels
# Moyenne des scores de dépendance par gène
essential_scores <- depmap %>%
  group_by(gene_name) %>%
  summarise(mean_dependency = mean(dependency, na.rm = TRUE)) %>%
  arrange(mean_dependency)

essential_scores <- depmap %>%  filter(dependency < -0.5)  %>%  arrange(dependency) 


library(tidyverse)
# Gènes essentiels colorectal (score < 0)
essential_genes <- essential_scores %>%
  filter(mean_dependency < -0) %>%
  pull(gene_name)


# Gènes essentiels colorectal (score <0)
essential_genes <- essential_scores %>%  pull(gene_name)


library(msigdbr)

load( "essential_genes_hepg2.rda")
load( "essential_scores_hepg2.rda")
load("pos.rda")
load("gene_list.rda")

# Récupérer les gènes oncogenic C6 humains
hallmark <- msigdbr(species = "Homo sapiens", category = "C5")


library(clusterProfiler)

# Préparer le vecteur nommé pour GSEA
gene_list <- pos$logFC  # ou autre colonne quantitative
names(gene_list) <- rownames(pos)
gene_list <- sort(gene_list, decreasing = TRUE)

# GSEA avec les Hallmarks
gsea_res <- GSEA(geneList = gene_list,
                 TERM2GENE = hallmark %>% dplyr::select(gs_name, gene_symbol),
                 pAdjustMethod = "BH",scoreType="pos",pvalueCutoff = 0.05,
                 verbose = FALSE)






dotplot(gsea_res, showCategory = 10)

gsea_df <- as.data.frame(gsea_res)

library(tidyr)

gsea_gene_map <- gsea_df %>%
  dplyr::select(ID = ID, Description = Description, core_enrichment) %>%
  separate_rows(core_enrichment, sep = "/") %>%
  dplyr::rename(Hallmark = Description, Gene = core_enrichment)


sl_gsea_annotated <- gsea_gene_map %>% filter(Gene %in% essential_genes)

library(igraph)
edges_gsea <- sl_gsea_annotated %>% dplyr::select(from = Gene, to = Hallmark)
g_gsea <- graph_from_data_frame(edges_gsea)

layout_opt <- layout_with_fr(g_gsea)
layout_opt <- layout_with_kk(g_gsea)  # Kamada-Kawai : plus aéré

plot(g_gsea,
     layout = layout_opt,
     vertex.label = V(g_gsea)$name,
     vertex.label.cex = 0.6,
     vertex.label.color = "black",
     vertex.color = ifelse(V(g_gsea)$name %in% essential_genes, "lightblue", "lightcoral"),
     vertex.frame.color = "gray",
     vertex.size = 10,
     vertex.shape = "circle",
     edge.arrow.size = 0.4,
     edge.color = "darkgray",
     main = "HB synthetic lethality –  GSEA")




###

library(ggraph)
library(tidygraph)

g_tbl <- as_tbl_graph(g_gsea)

ggraph(g_tbl, layout = "fr") +
  geom_edge_link(color = "gray70") +
  geom_node_point(aes(color = name %in% essential_genes), size = 4) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  theme_void()


save(essential_genes,file="essential_genes.rda")
save(essential_scores,file="essential_scores.rda")









library(dplyr)
library(tidyr)
library(purrr)

# Grouper les gènes par terme GSEA (Hallmark)
go_gene_groups <- sl_gsea_annotated %>%
  group_by(Hallmark) %>%
  summarise(Genes = list(unique(Gene)), .groups = "drop")

# Fonction pour générer les paires de gènes
generate_gene_pairs <- function(genes, term) {
  if (length(genes) < 2) return(NULL)
  pairs <- combn(genes, 2)
  data.frame(
    Gene1 = pairs[1, ],
    Gene2 = pairs[2, ],
    Term = term,
    stringsAsFactors = FALSE
  )
}

# Appliquer la fonction à chaque groupe
gene_pairs_df <- map2_dfr(go_gene_groups$Genes, go_gene_groups$Hallmark, generate_gene_pairs)


# Optionnel : retirer les doublons (Gene1-Gene2 vs Gene2-Gene1)
gene_pairs_df <- gene_pairs_df %>%
  mutate(pair_id = pmap_chr(list(Gene1, Gene2), ~ paste(sort(c(...)), collapse = "_"))) %>%
  distinct(pair_id, .keep_all = TRUE)


pos$gene<-row.names(pos)
## extract expression
expression_df <- pos %>%
  dplyr::select(gene, logFC) %>%
  dplyr::rename(Gene1 = gene, logFC_Gene1 = logFC)


# Renommer pour Gene2
expression_df2 <- pos %>%
  dplyr::select(gene, logFC) %>%
  dplyr::rename(Gene2 = gene, logFC_Gene2 = logFC)




# Calculer le score CRISPR moyen par gène
crispr_df <- depmap %>%
  group_by(gene_name) %>%
  summarise(mean_CRISPR = mean(dependency, na.rm = TRUE))

# Pour Gene1
crispr_gene1 <- crispr_df %>%
  dplyr::rename(Gene1 = gene_name, CRISPR_Gene1 = mean_CRISPR)

# Pour Gene2
crispr_gene2 <- crispr_df %>%
  dplyr::rename(Gene2 = gene_name, CRISPR_Gene2 = mean_CRISPR)

# Fusionner toutes les features
enriched_pairs <- gene_pairs_df %>%
  left_join(expression_df, by = "Gene1") %>%
  left_join(expression_df2, by = "Gene2") %>%
  left_join(crispr_gene1, by = "Gene1") %>%
  left_join(crispr_gene2, by = "Gene2")


# Définition empirique d'une bonne paire : logFC > 1 et CRISPR < -0.5 pour les deux gènes
enriched_pairs <- enriched_pairs %>%
  mutate(good_pair = ifelse(
    logFC_Gene1 > 1 & logFC_Gene2 > 1 & 
      CRISPR_Gene1 < -0.5 & CRISPR_Gene2 < -0.5,
    1, 0
  ))


library(randomForest)


library(dplyr)

# Transformer en tibble
enriched_pairs <- as_tibble(enriched_pairs)

enriched_pairs <- enriched_pairs %>%
  mutate(expression_score = (logFC_Gene1 + logFC_Gene2))


library(igraph)

# Fonction pour récupérer les termes GO d’un gène
get_go_terms <- function(gene, graph) {
  neighbors(graph, gene, mode = "out")$name
}

# Exemple : créer la colonne GO_similarity pour chaque paire
enriched_pairs$GO_similarity <- mapply(function(g1, g2) {
  terms1 <- get_go_terms(g1, g_gsea)
  terms2 <- get_go_terms(g2, g_gsea)
  inter <- length(intersect(terms1, terms2))
  union <- length(union(terms1, terms2))
  if (union == 0) return(0)
  return(inter / union) # Jaccard index
}, enriched_pairs$Gene1, enriched_pairs$Gene2)




# Degree
deg <- degree(g_gsea)
enriched_pairs$degree_Gene1 <- deg[enriched_pairs$Gene1]
enriched_pairs$degree_Gene2 <- deg[enriched_pairs$Gene2]



# Calcul des centralités
btw <- betweenness(g_gsea)
clo <- closeness(g_gsea)
eig <- eigen_centrality(g_gsea)$vector
pr  <- page.rank(g_gsea)$vector

# Nombre de voisins communs
enriched_pairs$common_neighbors <- mapply(function(g1, g2) {
  length(intersect(neighbors(g_gsea, g1), neighbors(g_gsea, g2)))
}, enriched_pairs$Gene1, enriched_pairs$Gene2)



# Centralités
btw <- betweenness(g_gsea)
clo <- closeness(g_gsea)
eig <- eigen_centrality(g_gsea)$vector
pr  <- page.rank(g_gsea)$vector

# Centralité combinée pour Gene1
enriched_pairs$centrality_G1 <- btw[enriched_pairs$Gene1] +
  clo[enriched_pairs$Gene1] +
  eig[enriched_pairs$Gene1] +
  pr[enriched_pairs$Gene1]

# Et Gene2 si besoin
enriched_pairs$centrality_G2 <- btw[enriched_pairs$Gene2] +
  clo[enriched_pairs$Gene2] +
  eig[enriched_pairs$Gene2] +
  pr[enriched_pairs$Gene2]



# Variante combinée : produit ou moyenne
enriched_pairs$centrality_sum <- enriched_pairs$centrality_G1 + enriched_pairs$centrality_G2
enriched_pairs$centrality_product <- enriched_pairs$centrality_G1 * enriched_pairs$centrality_G2

enriched_pairs$codep_score <- with(enriched_pairs, {
  # plus l'effet CRISPR est fort, plus l'inactivation est létale
  crispr_combo <- CRISPR_Gene1 * CRISPR_Gene2
  
  # centralité et interaction réseau
  net_factor <- log1p(centrality_product) + log1p(common_neighbors)
  
  # similarité fonctionnelle
  go_factor <- GO_similarity
  
  # optionnel : effet transcriptionnel
  perturbation <- abs(logFC_Gene1) + abs(logFC_Gene2)
  
  # score final (pondérations ajustables)
  (crispr_combo * 1.5) + (net_factor * 1.2) + (go_factor * 1.0) + (perturbation * 0.8)
})





features <- enriched_pairs[, c("logFC_Gene1", "logFC_Gene2",
                               "CRISPR_Gene1", "CRISPR_Gene2",
                               "GO_similarity", "degree_Gene1", "degree_Gene2",
                               "common_neighbors", "centrality_G1", "centrality_G2",
                               "centrality_sum", "centrality_product","expression_score")]
target <- enriched_pairs$codep_score  # ou codep_score selon l’objectif


library(randomForest)

set.seed(123)
rf_model <- randomForest(x = features, y = target, ntree = 500, importance = TRUE)

# Score prédits
enriched_pairs$AI_score <- predict(rf_model, features)

top_ranked_pairs <- enriched_pairs[order(-enriched_pairs$AI_score), ]
head(top_ranked_pairs[, c("Gene1", "Gene2", "GO_similarity", "AI_score")], 30)


varImpPlot(rf_model)          # Visualiser les features les plus influents
importance(rf_model)          # Valeurs numériques d’importance


library(ggplot2)
library(dplyr)

# Créer une colonne "Pair" pour fusionner les noms des gènes
top_ranked_pairs <- top_ranked_pairs %>%
  mutate(Pair = paste(Gene1, Gene2, sep = "–"))



top_ranked_pairs <- top_ranked_pairs %>%
  mutate(
    Gene1 = as.character(Gene1),
    Gene2 = as.character(Gene2),
    GO_BP = as.character(GO_similarity),
    AI_score = as.numeric(AI_score),
    Pair = paste(Gene1, Gene2, sep = "–")
  )


# Représentation ggplot
top_ranked_pairs%>%head(.,n=35)%>%
ggplot(., aes(x = reorder(Pair, AI_score), y = AI_score, fill = Term)) +
  geom_col() +
  coord_flip() +
  labs(title = "",
       x = "Pairs of genes",
       y = "RF_importance_score",
       fill = "Term") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")

plot(rf_model)

sl_go_annotated%>%filter(Gene=="CCNB1")


save(enriched_pairs,file="enriched_pairs.rda")

write.csv(top_ranked_pairs,file="topPairs.csv")


## export cytoscape
write_graph(g_gsea, file = "g_go_network.graphml", format = "graphml")

## rf interpretation

preds <- predict(rf_model)
plot(target, preds, xlab = "real values", ylab = "predictions", main = "Prediction vs Real values")
abline(0, 1, col = "red", lwd = 2)  # ligne parfaite

r_squared <- 1 - (sum((target - preds)^2) / sum((target - mean(target))^2))



library(caret)
ctrl <- trainControl(method = "cv", number = 5)
rf_cv <- train(target ~ ., data = data.frame(features, target),
               method = "rf", trControl = ctrl)
rf_cv



save(sl_gsea_annotated,file="sl_gsea_annotated.rda")
write.table(sl_gsea_annotated,file="sl_gsea_annotated.tsv",row.names=F)

sl_gsea_annotated%>%distinct(Gene)%>%pull(Gene)->input

length(input)

load("quantile.rda")

pheno<-read.csv("pheno.csv",h=T,row.names=1)

data<-norm_edata
all(colnames(data)==row.names(pheno))

df2<-data[row.names(data)%in%input,]

pheno%>%select(group,batch)->annot2
library(transpipe15)
pheatmap(df2, scale = "row", color = colorRampPalette(c("navy", 
            "white", "firebrick3"))(50), annotation = annot2, 
            fontsize = 12, cutree_rows = 1, cutree_col = 1, 
            clustering_method = "ward.D2", clustering_distance_cols = "euclidean", 
            clustering_distance_rows = "euclidean", show_colnames = F, cluster_row=F,
            show_rownames = T, annotation_names_col = T, annotation_names_row = T)


depmap%>%filter(gene_name%in%input)->df


library(tidyr)
library(dplyr)

# Exemple : df est ton tableau initial
df_wide <- df %>%
  select(gene_name, cell_line, dependency) %>%
  pivot_wider(names_from = cell_line, values_from = dependency)

df<-as.data.frame(df_wide)
row.names(df)<-df$gene_name
df$gene_name<-NULL
bestheat(df,annot,scale="row",font=10,rownames=T)


annot<-as.data.frame(colnames(df))
colnames(annot)<-"id"
row.names(annot)<-annot$id
pheatmap(df, scale = "none", color = colorRampPalette(c("firebrick3","white"))(50), annotation = annot, 
            fontsize = 12, cutree_rows = 1, cutree_col = 1, 
            clustering_method = "ward.D2", clustering_distance_cols = "euclidean", 
            clustering_distance_rows = "euclidean", show_colnames = F, cluster_row=F,
            show_rownames = T, annotation_names_col = T, annotation_names_row = T)


## centrality
centrality_df <- data.frame(
  Node = names(btw),
  Betweenness = btw,
  Closeness = clo[names(btw)],
  Eigenvector = eig[names(btw)],
  PageRank = pr[names(btw)],
  stringsAsFactors = FALSE
)


centrality_df %>% arrange(desc(PageRank)) ->centrality

library(stringr)
centrality %>%  filter(Closeness!="NaN")%>% arrange(desc(Eigenvector))->centrality

ggplot(centrality, aes(x = Node, y = Eigenvector)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +  # ← rend le barplot vertical (nœuds en Y)
  labs(title = "Eigenvector centrality of nodes",
       x = "Nodes",
       y = "Score of centrality") +
  theme_minimal(base_size = 12)



save(centrality,file="centrality.rda")


load("centrality.rda")
library(dplyr)
centrality%>%arrange(desc(Node))->centrality


library(dplyr)

centrality_ranked <- centrality %>%
  arrange(desc(Node)) %>%  # Tri décroissant sur Eigenvector
  mutate(Index = row_number()) %>%
  select(Index, everything())     # Place l'index en première colonne


library(ggplot2)




centrality_ranked <- centrality_ranked %>%
  mutate(Index = as.numeric(Index),
         Node = reorder(Node, Index))

ggplot(centrality_ranked, aes(x = Node, y = Eigenvector)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Eigenvector centrality",
       x = "Nodes",
       y = "Score of centrality") +
  theme_minimal(base_size = 12)

write.csv(centrality,file="gene_centrality_42.csv",row.names=F)
