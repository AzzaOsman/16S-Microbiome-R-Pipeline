River_16S_Analysis.R (with Unfiltered + Filtered branches)
Step 19
###############################################################################
# River Microbiome Analysis (16S rRNA)
# Author: Azzah
# Updated: 2025-11-13
# Description: Full DADA2 + Phyloseq pipeline with unfiltered alpha diversity
#              and filtered beta diversity / relative abundance / heatmap
###############################################################################
# 19. BETA DIVERSITY (Filtered)
###############################################################################

bray_dist <- phyloseq::distance(ps_filtered, method = "bray")
metadata_df <- data.frame(sample_data(ps_filtered))

# PERMANOVA
adonis_result <- adonis2(bray_dist ~ River, data = metadata_df)
adonis_df <- as.data.frame(adonis_result)

# Pairwise PERMANOVA
pairwise_result <- pairwise.adonis(bray_dist, factors = metadata_df$River, sim.method = "bray", p.adjust.m = "bonferroni")
pairwise_df <- as.data.frame(pairwise_result)

# Save PERMANOVA results to Excel
write.xlsx(list(
  PERMANOVA = adonis_df,
  Pairwise_PERMANOVA = pairwise_df
), file = file.path(result_path, "BetaDiversity_PERMANOVA.xlsx"))

library(ggrepel)

# PCoA with better labels
p_pcoa <- plot_ordination(ps_filtered, ordination_pcoa, color = "River") +
  geom_point(size = 4) +
  geom_text_repel(aes(label = sample_names(ps_filtered)), size = 3, max.overlaps = 50) +
  theme_bw() +
  labs(title = "PCoA (Bray-Curtis)", x = "PCoA1", y = "PCoA2")

ggsave(file.path(result_path, "PCoA_BrayCurtis_Labelled_Improved.png"),
       plot = p_pcoa, width = 10, height = 8, dpi = 300)

# NMDS with better labels
p_nmds <- plot_ordination(ps_filtered, ordination_nmds, color = "River") +
  geom_point(size = 4) +
  geom_text_repel(aes(label = sample_names(ps_filtered)), size = 3, max.overlaps = 50) +
  theme_bw() +
  labs(title = "NMDS (Bray-Curtis)")

ggsave(file.path(result_path, "NMDS_BrayCurtis_Labelled_Improved.png"),
       plot = p_nmds, width = 10, height = 8, dpi = 300)
