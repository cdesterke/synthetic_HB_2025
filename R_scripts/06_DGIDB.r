library(dplyr)
# load drug gene interaction database
interactions <- read.delim("interactions.tsv")
categories <- read.delim("categories.tsv")

# load HB essential genes
genes_of_interest<-c("FABP4","PEG10","PLCB4","SEMA7A","HDAC11","ZNF233")

interactions_clean <- interactions %>%
  filter(!is.na(interaction_type))

## select repressive drugs
repressive_types <- c("inhibitor", "antagonist", "blocker", "suppressor",
                      "negative modulator", "inverse agonist", "cleavage",
                      "antisense oligonucleotide")

repressive_interactions <- interactions %>%
  filter(interaction_type %in% repressive_types,
         gene_name %in% genes_of_interest) %>%
  select(gene_name, drug_name, interaction_type, interaction_score,approved)

write.csv(repressive_interactions,file="repressive_interactions.csv",row.names=F)


## barplots
library(dplyr)
library(ggplot2)

# Regrouper par gène et déterminer si au moins un médicament est approuvé
gene_drug_summary <- repressive_interactions %>%
  group_by(gene_name, approved) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(approval_status = if_else(approved == TRUE, "Approved", "Not approved")) %>%
  arrange(desc(n))


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


repressive_interactions %>%
  mutate(interaction_score = as.numeric(as.character(interaction_score)))


library(ggplot2)
library(dplyr)

# prepare data
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


write.csv(repressive_summary,file="NEW_actionable.csv",row.names=F)


interactions <- interactions %>%
  mutate(interaction_score = as.numeric(interaction_score))

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




## network
library(igraph)
library(ggraph)
library(tidyverse)




### type of interactions




edges <- repressive_interactions %>%
  mutate(interaction_score = as.numeric(interaction_score)) %>%
  select(from = drug_name, to = gene_name,
         type = interaction_type, weight = interaction_score)

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
  labs(
    title = "Network of repressive drug interactions (DGIdb)",
    subtitle = "Type of interactions (color) – genes of interest ↔ Drugs",
    edge_color = "Type of interaction",
    edge_width = "Score of interaction"
  ) +
  theme_void()







## network with drug approbation
# Edges: from drug to gene
library(tidyverse)
library(igraph)
library(ggraph)

# ---- 1. Edges ----
edges <- repressive_interactions %>%
  transmute(
    from   = drug_name,
    to     = gene_name,
    type   = interaction_type,
    weight = as.numeric(interaction_score)
  )

# ---- 2. Nodes ----
nodes <- tibble(name = unique(c(edges$from, edges$to))) %>%
  mutate(
    node_type = if_else(name %in% edges$from, "Drug", "Gene"),
    approved  = case_when(
      node_type == "Gene" ~ "Not applicable",
      name %in% repressive_interactions$drug_name[repressive_interactions$approved == TRUE] ~ "Approved drug",
      TRUE ~ "Not approved drug"
    ),
    shape = case_when(
      node_type == "Gene" ~ "gene_shape",              # forme spécifique pour gènes
      approved == "Approved drug" ~ "approved_shape",  # cercle
      approved == "Not approved drug" ~ "not_approved_shape"  # carré
    )
  )

# ---- 3. Graph ----
g <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)

# ---- 4. Plot ----
ggraph(g, layout = "fr") +
  geom_edge_link(aes(color = type, width = weight), alpha = 1) +
  geom_node_point(aes(color = node_type, shape = shape), size = 5) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +

  # couleurs des interactions
  scale_edge_color_manual(values = c(
    "inhibitor" = "green",
    "cleavage"  = "blue"
  )) +

  # couleurs des nœuds
  scale_color_manual(values = c(
    "Drug" = "#E64B35",
    "Gene" = "#4DBBD5"
  )) +

  # formes explicites
  scale_shape_manual(values = c(
    "gene_shape"          = 17,  # triangle pour gènes
    "approved_shape"      = 16,  # cercle pour drugs approuvés
    "not_approved_shape"  = 15   # carré pour drugs non approuvés
  )) +

  scale_edge_width(range = c(0.2, 2)) +

  labs(
    title = "Repressive Drug–Gene Interaction Network",
    subtitle = "Node shape indicates drug approval status",
    edge_color = "Interaction Type",
    color = "Node Type",
    shape = "Approval Status"
  ) +
  theme_void()






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



repressive_interactions%>%filter(approved=="TRUE")->approved
