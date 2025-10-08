library(dplyr)
# Exemple de croisement
interactions <- read.delim("interactions.tsv")
categories <- read.delim("categories.tsv")

net<-read.csv("gene_centrality_42.csv",h=T)
genes_of_interest<-net$Node

interactions_clean <- interactions %>%
  filter(!is.na(interaction_type))

repressive_types <- c("inhibitor", "antagonist", "blocker", "suppressor",
                      "negative modulator", "inverse agonist", "cleavage",
                      "antisense oligonucleotide")

repressive_interactions <- interactions %>%
  filter(interaction_type %in% repressive_types,
         gene_name %in% genes_of_interest) %>%
  select(gene_name, drug_name, interaction_type, interaction_score,approved)

write.csv(repressive_interactions,file="repressive_interactions.csv",row.names=F)


library(dplyr)
library(ggplot2)

# Regrouper par gène et déterminer si au moins un médicament est approuvé
gene_drug_summary <- repressive_interactions %>%
  group_by(gene_name, approved) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(approval_status = if_else(approved == TRUE, "Approved", "Not approved")) %>%
  arrange(desc(n))



# Barplot
ggplot(gene_drug_summary, aes(x = reorder(gene_name, n), y = n, fill = approval_status)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Approved" = "forestgreen", "Not approved" = "firebrick")) +
  labs(title = "Number of Repressive Drugs per Gene",
       subtitle = "Colored by Approval Status",
       x = "Genes",
       y = "Number of Drugs",
       fill = "Approval Status") +
  theme_minimal(base_size = 12)
















repressive_interactions %>%
  count(gene_name, interaction_type) %>% arrange(desc(n))


 <- repressive_interactions %>%
  mutate(interaction_score = as.numeric(as.character(interaction_score)))


library(ggplot2)
library(dplyr)

# Données préparées
repressive_summary <- repressive_interactions %>%
  count(gene_name, interaction_type) %>%
  arrange(desc(n))

# Barplot
ggplot(repressive_summary, aes(x = reorder(gene_name, n), y = n, fill = interaction_type)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Repressive interactions by gene (DGIdb)",
       x = "Genes",
       y = "Numbers of interactions",
       fill = "Type of interactions") +
  theme_minimal(base_size = 12) +
  scale_fill_manual(values = c("inhibitor" = "green", "cleavage" = "blue")) 


write.csv(repressive_summary,file="eleven_actionable.csv",row.names=F)



repressive_by_drug <- interactions %>%
  filter(interaction_type %in% repressive_types,
         gene_name %in% genes_of_interest) %>%
  group_by(drug_name) %>%
  summarise(
    n_genes_targeted = n_distinct(gene_name),
    genes = paste(unique(gene_name), collapse = ", "),
    avg_score = round(mean(interaction_score, na.rm = TRUE), 4)
  ) %>%
  arrange(desc(n_genes_targeted))



library(igraph)
library(ggraph)
library(tidyverse)




### type of interactions

# Créer les arêtes avec type et score
edges <- repressive_interactions %>%
  select(from = drug_name, to = gene_name, type = interaction_type, weight = interaction_score)

# Créer les nœuds avec type (gène ou médicament)
nodes <- tibble(name = unique(c(edges$from, edges$to))) %>%
  mutate(node_type = if_else(name %in% edges$from, "Drug", "Gene"))


g <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)


ggraph(g, layout = "fr") +
  geom_edge_link(aes(color = type, width = weight), alpha = 1) +
  geom_node_point(aes(color = node_type), size = 5) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_edge_color_manual(values = c(
    "inhibitor" = "green",
    "antagonist" = "#FC8D59",
    "blocker" = "#91BFDB",
    "suppressor" = "#4575B4",
    "negative modulator" = "#313695",
    "inverse agonist" = "#74ADD1",
    "cleavage" = "blue",
    "antisense oligonucleotide" = "#E0F3F8"
  )) +
  scale_color_manual(values = c("Drug" = "#E64B35", "Gene" = "#4DBBD5")) +
  scale_edge_width(range = c(0.2, 2)) +
  labs(title = "Network of repressive drug interactions (DGIdb)",
       subtitle = "Type of interactions (color) – genes of interest ↔ Drugs",
       edge_color = "Type of interaction",
       edge_width = "Score of interaction") +
  theme_void()






## net approbation
# Edges: from drug to gene
edges <- repressive_interactions %>%
  select(from = drug_name, to = gene_name, type = interaction_type, weight = interaction_score)

# Nodes: combine drug and gene names
nodes <- tibble(name = unique(c(edges$from, edges$to))) %>%
  mutate(
    node_type = if_else(name %in% edges$from, "Drug", "Gene"),
    approved = if_else(node_type == "Drug" & name %in% repressive_interactions$drug_name[repressive_interactions$approved == TRUE], "Approved", "Not approved"),
    shape = case_when(
      node_type == "Gene" ~ "circle",
      approved == "Approved" ~ "circle",
      approved == "Not approved" ~ "square",
      TRUE ~ "circle"
    )
  )

g <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)


ggraph(g, layout = "fr") +
  geom_edge_link(aes(color = type, width = weight), alpha = 1) +
  geom_node_point(aes(color = node_type, shape = shape), size = 5) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_edge_color_manual(values = c(
    "inhibitor" = "green",     # forest green
    "cleavage" = "blue"       # royal blue
  )) +
  scale_color_manual(values = c("Drug" = "#E64B35", "Gene" = "#4DBBD5")) +
  scale_shape_manual(values = c("circle" = 16, "square" = 15)) +
  scale_edge_width(range = c(0.2, 2)) +
  labs(title = "Repressive Drug–Gene Interaction Network",
       subtitle = "Node shape indicates drug approval status",
       edge_color = "Interaction Type",
       shape = "Approval Status",
       color = "Node Type") +
  theme_void()



repressive_interactions %>%filter(gene_name=="ATP1A2")


drug_type_summary <- repressive_interactions %>%
  group_by(gene_name, interaction_type) %>%
  summarise(n_drugs = n_distinct(drug_name), .groups = "drop")


ggplot(drug_type_summary, aes(x = reorder(gene_name, n_drugs), y = n_drugs, fill = interaction_type)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +
  scale_fill_manual(values = c("inhibitor" = "green", "cleavage" = "blue")) +
  labs(title = "Unique Repressive Drugs per Gene",
       subtitle = "Stacked by Interaction Type",
       x = "Genes",
       y = "Number of Unique Drugs",
       fill = "Interaction Type") +
  theme_minimal(base_size = 12)




