# data of dependency from depmap R-package

load("depmap_scores.rda")
# metadata from cell lines

load("depmap_samples.rda")


# Extraction of hepatoblastoma specific data
library(dplyr)
library(stringr)
lines <-sample_info %>% filter(str_detect(Cellosaurus_NCIt_disease, regex("Hepatoblastoma", ignore_case = TRUE)))
lines<-as.data.frame(lines)
write.csv(sample_info,file="metadata.csv",row.names=F)

lines%>%  pull(depmap_id) %>% unique()
id<-c("ACH-000739","ACH-000671")

# filtration of CRISPR scores for HB cell lines
depmap <- depmap_scores %>%
  filter(depmap_id %in% id)

## identification of essential genes
essential_scores <- depmap %>%
  group_by(gene_name) %>%
  summarise(mean_dependency = mean(dependency, na.rm = TRUE)) %>%
  arrange(mean_dependency)

essential_genes <- essential_scores %>%
  filter(mean_dependency < 0) %>%
  pull(gene_name)


## GSEA GO
library(msigdbr)
## load limma up regulated genes
load("pos.rda")


# GSEA GO genesets
hallmark <- msigdbr(species = "Homo sapiens", category = "C5")


library(clusterProfiler)

# prepare data GSEA
gene_list <- pos$logFC  # quantitative column
names(gene_list) <- rownames(pos)
gene_list <- sort(gene_list, decreasing = TRUE)

# GSEA 
gsea_res <- GSEA(geneList = gene_list,
                 TERM2GENE = hallmark %>% dplyr::select(gs_name, gene_symbol),
                 pAdjustMethod = "BH",scoreType="pos",pvalueCutoff = 0.05,
                 verbose = FALSE)

dotplot(gsea_res, showCategory = 10)


## extract data from enrichment
gsea_df <- as.data.frame(gsea_res)
library(tidyr)
gsea_gene_map <- gsea_df %>%
  dplyr::select(ID = ID, Description = Description, core_enrichment) %>%
  separate_rows(core_enrichment, sep = "/") %>%
  dplyr::rename(Hallmark = Description, Gene = core_enrichment)

## filter the enrichment with essential genes
sl_gsea_annotated <- gsea_gene_map %>% filter(Gene %in% essential_genes)

library(igraph)
edges_gsea <- sl_gsea_annotated %>% dplyr::select(from = Gene, to = Hallmark)
g_gsea <- graph_from_data_frame(edges_gsea)


### tidy network

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

# group genes by GSEA terms
go_gene_groups <- sl_gsea_annotated %>%
  group_by(Hallmark) %>%
  summarise(Genes = list(unique(Gene)), .groups = "drop")

# generation of pairs of genes
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

# applied generate_gene_pairs() to each groups
gene_pairs_df <- map2_dfr(go_gene_groups$Genes, go_gene_groups$Hallmark, generate_gene_pairs)


# remove doublons (Gene1-Gene2 vs Gene2-Gene1)
gene_pairs_df <- gene_pairs_df %>%
  mutate(pair_id = pmap_chr(list(Gene1, Gene2), ~ paste(sort(c(...)), collapse = "_"))) %>%
  distinct(pair_id, .keep_all = TRUE)


## extract expression
pos$gene<-row.names(pos)
expression_df <- pos %>%
  dplyr::select(gene, logFC) %>%
  dplyr::rename(Gene1 = gene, logFC_Gene1 = logFC)


# rename Gene2
expression_df2 <- pos %>%
  dplyr::select(gene, logFC) %>%
  dplyr::rename(Gene2 = gene, logFC_Gene2 = logFC)

#Compute mean CRISPR score by gene
crispr_df <- depmap %>%
  group_by(gene_name) %>%
  summarise(mean_CRISPR = mean(dependency, na.rm = TRUE))

# for Gene1
crispr_gene1 <- crispr_df %>%
  dplyr::rename(Gene1 = gene_name, CRISPR_Gene1 = mean_CRISPR)

# for Gene2
crispr_gene2 <- crispr_df %>%
  dplyr::rename(Gene2 = gene_name, CRISPR_Gene2 = mean_CRISPR)

# aggegates
enriched_pairs <- gene_pairs_df %>%
  left_join(expression_df, by = "Gene1") %>%
  left_join(expression_df2, by = "Gene2") %>%
  left_join(crispr_gene1, by = "Gene1") %>%
  left_join(crispr_gene2, by = "Gene2")


# definition of filters
enriched_pairs <- enriched_pairs %>%
  mutate(good_pair = ifelse(
    logFC_Gene1 > 1 & logFC_Gene2 > 1 & 
      CRISPR_Gene1 < -0.5 & CRISPR_Gene2 < -0.5,
    1, 0
  ))


library(dplyr)

# transform in tibble
enriched_pairs <- as_tibble(enriched_pairs)

enriched_pairs <- enriched_pairs %>%
  mutate(expression_score = (logFC_Gene1 + logFC_Gene2))


library(igraph)

# Go terms by gene
get_go_terms <- function(gene, graph) {
  neighbors(graph, gene, mode = "out")$name
}

# Go similarity between 2 genes
enriched_pairs$GO_similarity <- mapply(function(g1, g2) {
  terms1 <- get_go_terms(g1, g_gsea)
  terms2 <- get_go_terms(g2, g_gsea)
  inter <- length(intersect(terms1, terms2))
  union <- length(union(terms1, terms2))
  if (union == 0) return(0)
  return(inter / union) # Jaccard index
}, enriched_pairs$Gene1, enriched_pairs$Gene2)

# Degree of network
deg <- degree(g_gsea)
enriched_pairs$degree_Gene1 <- deg[enriched_pairs$Gene1]
enriched_pairs$degree_Gene2 <- deg[enriched_pairs$Gene2]

# compute centrality
btw <- betweenness(g_gsea)
clo <- closeness(g_gsea)
eig <- eigen_centrality(g_gsea)$vector
pr  <- page.rank(g_gsea)$vector

# number of neighbors
enriched_pairs$common_neighbors <- mapply(function(g1, g2) {
  length(intersect(neighbors(g_gsea, g1), neighbors(g_gsea, g2)))
}, enriched_pairs$Gene1, enriched_pairs$Gene2)



# Combined centrality for Gene 1
enriched_pairs$centrality_G1 <- btw[enriched_pairs$Gene1] +
  clo[enriched_pairs$Gene1] +
  eig[enriched_pairs$Gene1] +
  pr[enriched_pairs$Gene1]

# and for Gene 2
enriched_pairs$centrality_G2 <- btw[enriched_pairs$Gene2] +
  clo[enriched_pairs$Gene2] +
  eig[enriched_pairs$Gene2] +
  pr[enriched_pairs$Gene2]

# combinations of centrality
enriched_pairs$centrality_sum <- enriched_pairs$centrality_G1 + enriched_pairs$centrality_G2
enriched_pairs$centrality_product <- enriched_pairs$centrality_G1 * enriched_pairs$centrality_G2

# codependency scores
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

## randomForest to predict codep_scores with variables of the map
library(randomForest)

set.seed(123)
rf_model <- randomForest(x = features, y = target, ntree = 500, importance = TRUE)

# predicted scores
enriched_pairs$AI_score <- predict(rf_model, features)

# ranking on predicted score
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


# ggplot2  representation of the ranked pairs
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


save(enriched_pairs,file="enriched_pairs.rda")

write.csv(top_ranked_pairs,file="topPairs.csv")


## export cytoscape
write_graph(g_gsea, file = "g_go_network.graphml", format = "graphml")

## rf interpretation

preds <- predict(rf_model)
plot(target, preds, xlab = "real values", ylab = "predictions", main = "Prediction vs Real values")
abline(0, 1, col = "red", lwd = 2)  # ligne parfaite

r_squared <- 1 - (sum((target - preds)^2) / sum((target - mean(target))^2))





save(sl_gsea_annotated,file="sl_gsea_annotated.rda")
write.table(sl_gsea_annotated,file="sl_gsea_annotated.tsv",row.names=F)


## heatmap visualization
sl_gsea_annotated%>%distinct(Gene)%>%pull(Gene)->input

length(input)

## normalized expression data
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

# table wider
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


## barplot centrality
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

