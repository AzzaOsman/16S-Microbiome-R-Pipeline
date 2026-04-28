River_16S_Analysis.R (with Unfiltered + Filtered branches)
Step 18
###############################################################################
# River Microbiome Analysis (16S rRNA)
# Author: Azza
# Updated: 2025-11-13
# Description: Full DADA2 + Phyloseq pipeline with unfiltered alpha diversity
#              and filtered beta diversity / relative abundance / heatmap
###############################################################################
# Custom colours for River annotation
river_colors <- c("Jelau" = "#E41A1C",       # Red
                  "Lubok Paoh" = "#377EB8",  # Blue
                  "Ng Linsum" = "#4DAF4A",   # Green
                  "Ng Spak" = "#984EA3",     # Purple
                  "Ng Tiga" = "#FF7F00")     # Orange

library(pheatmap)

#Phylum-level heatmap
# Phylum-level data
ps_phylum <- tax_glom(ps_filtered, taxrank = "Phylum")
ps_phylum_rel <- transform_sample_counts(ps_phylum, function(x) x / sum(x))

# Extract abundance matrix
heatmap_mat_phylum <- as.data.frame(otu_table(ps_phylum_rel))
if (!taxa_are_rows(ps_phylum_rel)) heatmap_mat_phylum <- t(heatmap_mat_phylum)

# Replace row names with Phylum names
tax_names_phylum <- tax_table(ps_phylum)[, "Phylum"] %>% as.character()
rownames(heatmap_mat_phylum) <- tax_names_phylum

# Annotation for samples
sample_info <- data.frame(River = sample_data(ps_phylum)$River)
rownames(sample_info) <- sample_names(ps_phylum)

# Annotation colours
ann_colors <- list(River = river_colors)

# Heatmap
pheatmap(heatmap_mat_phylum,
         cluster_rows = TRUE, cluster_cols = TRUE,
         main = "Phylum-Level Heatmap",
         annotation_col = sample_info,
         annotation_colors = ann_colors,
         angle_col = 45,
         fontsize_row = 10, fontsize_col = 12, fontsize = 14,
         filename = file.path(path, "Result/Phylum_Heatmap_Annotated.png"),
         width = 12, height = 10)

#Family-level Heatmap (Top 30)
# Family-level data
ps_family <- tax_glom(ps_filtered, taxrank = "Family")
ps_family_rel <- transform_sample_counts(ps_family, function(x) x / sum(x))

heatmap_mat_family <- as.data.frame(otu_table(ps_family_rel))
if (!taxa_are_rows(ps_family_rel)) heatmap_mat_family <- t(heatmap_mat_family)

tax_names_family <- tax_table(ps_family)[, "Family"] %>% as.character()
rownames(heatmap_mat_family) <- tax_names_family

# Top 30 families
top_families <- rowSums(heatmap_mat_family) %>%
  sort(decreasing = TRUE) %>%
  head(30) %>%
  names()

heatmap_mat_family <- as.matrix(heatmap_mat_family[top_families, ])

# Annotation
sample_info <- data.frame(River = sample_data(ps_family)$River)
rownames(sample_info) <- sample_names(ps_family)

pheatmap(heatmap_mat_family,
         cluster_rows = TRUE, cluster_cols = TRUE,
         main = "Family-Level Heatmap (Top 30)",
         annotation_col = sample_info,
         annotation_colors = ann_colors,
         angle_col = 45,
         fontsize_row = 10, fontsize_col = 12, fontsize = 14,
         filename = file.path(path, "Result/Family_Heatmap_Top30_Annotated.png"),
         width = 12, height = 10)

#Genus-level heatmap (Top 30)
# Genus-level data
ps_genus <- tax_glom(ps_filtered, taxrank = "Genus")
ps_genus_rel <- transform_sample_counts(ps_genus, function(x) x / sum(x))

heatmap_mat_genus <- as.data.frame(otu_table(ps_genus_rel))
if (!taxa_are_rows(ps_genus_rel)) heatmap_mat_genus <- t(heatmap_mat_genus)

tax_names_genus <- tax_table(ps_genus)[, "Genus"] %>% as.character()
rownames(heatmap_mat_genus) <- tax_names_genus

# Top 30 genera
top_genera <- rowSums(heatmap_mat_genus) %>%
  sort(decreasing = TRUE) %>%
  head(30) %>%
  names()

heatmap_mat_genus <- as.matrix(heatmap_mat_genus[top_genera, ])

# Annotation
sample_info <- data.frame(River = sample_data(ps_genus)$River)
rownames(sample_info) <- sample_names(ps_genus)

pheatmap(heatmap_mat_genus,
         cluster_rows = TRUE, cluster_cols = TRUE,
         main = "Genus-Level Heatmap (Top 30)",
         annotation_col = sample_info,
         annotation_colors = ann_colors,
         angle_col = 45,
         fontsize_row = 10, fontsize_col = 12, fontsize = 14,
         filename = file.path(path, "Result/Genus_Heatmap_Top30_Annotated.png"),
         width = 12, height = 10)

#Heatmap with colour gradient blue-white-red for abundance
#Define River Colour and Gradient
# Custom colours for River annotation
river_colors <- c("Jelau" = "#E41A1C",       # Red
                  "Lubok Paoh" = "#377EB8",  # Blue
                  "Ng Linsum" = "#4DAF4A",   # Green
                  "Ng Spak" = "#984EA3",     # Purple
                  "Ng Tiga" = "#FF7F00")     # Orange

# Annotation colour list
ann_colors <- list(River = river_colors)

# Colour gradient for abundance
heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)

#Phylum Level Heatmap
ps_phylum <- tax_glom(ps_filtered, taxrank = "Phylum")
ps_phylum_rel <- transform_sample_counts(ps_phylum, function(x) x / sum(x))

heatmap_mat_phylum <- as.data.frame(otu_table(ps_phylum_rel))
if (!taxa_are_rows(ps_phylum_rel)) heatmap_mat_phylum <- t(heatmap_mat_phylum)

tax_names_phylum <- tax_table(ps_phylum)[, "Phylum"] %>% as.character()
rownames(heatmap_mat_phylum) <- tax_names_phylum

sample_info <- data.frame(River = sample_data(ps_phylum)$River)
rownames(sample_info) <- sample_names(ps_phylum)

pheatmap(heatmap_mat_phylum,
         cluster_rows = TRUE, cluster_cols = TRUE,
         main = "Phylum-Level Heatmap",
         annotation_col = sample_info,
         annotation_colors = ann_colors,
         color = heatmap_colors,
         angle_col = 45,
         fontsize_row = 10, fontsize_col = 12, fontsize = 14,
         filename = file.path(path, "Result/Phylum_Heatmap_Annotated.png"),
         width = 12, height = 10)

#Family-level heatmap (Top 30)

ps_family <- tax_glom(ps_filtered, taxrank = "Family")
ps_family_rel <- transform_sample_counts(ps_family, function(x) x / sum(x))

heatmap_mat_family <- as.data.frame(otu_table(ps_family_rel))
if (!taxa_are_rows(ps_family_rel)) heatmap_mat_family <- t(heatmap_mat_family)

tax_names_family <- tax_table(ps_family)[, "Family"] %>% as.character()
rownames(heatmap_mat_family) <- tax_names_family

top_families <- rowSums(heatmap_mat_family) %>%
  sort(decreasing = TRUE) %>%
  head(30) %>%
  names()

heatmap_mat_family <- as.matrix(heatmap_mat_family[top_families, ])

sample_info <- data.frame(River = sample_data(ps_family)$River)
rownames(sample_info) <- sample_names(ps_family)

pheatmap(heatmap_mat_family,
         cluster_rows = TRUE, cluster_cols = TRUE,
         main = "Family-Level Heatmap (Top 30)",
         annotation_col = sample_info,
         annotation_colors = ann_colors,
         color = heatmap_colors,
         angle_col = 45,
         fontsize_row = 10, fontsize_col = 12, fontsize = 14,
         filename = file.path(path, "Result/Family_Heatmap_Top30_Annotated.png"),
         width = 12, height = 10)

#Genus-level Heatmap (Top 30)
ps_genus <- tax_glom(ps_filtered, taxrank = "Genus")
ps_genus_rel <- transform_sample_counts(ps_genus, function(x) x / sum(x))

heatmap_mat_genus <- as.data.frame(otu_table(ps_genus_rel))
if (!taxa_are_rows(ps_genus_rel)) heatmap_mat_genus <- t(heatmap_mat_genus)

tax_names_genus <- tax_table(ps_genus)[, "Genus"] %>% as.character()
rownames(heatmap_mat_genus) <- tax_names_genus

top_genera <- rowSums(heatmap_mat_genus) %>%
  sort(decreasing = TRUE) %>%
  head(30) %>%
  names()

heatmap_mat_genus <- as.matrix(heatmap_mat_genus[top_genera, ])

sample_info <- data.frame(River = sample_data(ps_genus)$River)
rownames(sample_info) <- sample_names(ps_genus)

pheatmap(heatmap_mat_genus,
         cluster_rows = TRUE, cluster_cols = TRUE,
         main = "Genus-Level Heatmap (Top 30)",
         annotation_col = sample_info,
         annotation_colors = ann_colors,
         color = heatmap_colors,
         angle_col = 45,
         fontsize_row = 10, fontsize_col = 12, fontsize = 14,
         filename = file.path(path, "Result/Genus_Heatmap_Top30_Annotated.png"),
         width = 12, height = 10)
