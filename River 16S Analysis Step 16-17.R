River_16S_Analysis.R (with Unfiltered + Filtered branches)
Step 16-17
###############################################################################
# River Microbiome Analysis (16S rRNA)
# Author: Azza
# Updated: 2025-11-13
# Description: Full DADA2 + Phyloseq pipeline with unfiltered alpha diversity
#              and filtered beta diversity / relative abundance / heatmap
###############################################################################
# 16. CREATE FILTERED PHYLOSEQ OBJECT (for Beta, Relative Abundance, Heatmap)
###############################################################################
table(rowSums(otu_table(ps_unfiltered) > 0))
ps_filtered <- ps_unfiltered  # skip filtering if too strict

# Moderate filtering: keep taxa present in ≥2 samples and total count ≥10
ps_filtered <- filter_taxa(ps_filtered, function(x) sum(x > 0) >= 2, prune = TRUE)

ps_filtered <- prune_taxa(taxa_sums(ps_filtered) >= 10, ps_filtered)

saveRDS(ps_filtered, "C:/Users/erino/Desktop/13-11-2025 River R/Result/ps_filtered_for_beta_relative_heatmap.rds")

###############################################################################
# 17. RELATIVE ABUNDANCE PLOTS (Phylum, Family, Genus)
###############################################################################
# --- Phylum Level ---
if (!dir.exists(file.path(path, "Result"))) dir.create(file.path(path, "Result"), recursive = TRUE)

# --- Phylum Level ---
# Prepare Phylum-level data
ps_phylum <- tax_glom(ps_filtered, taxrank = "Phylum")
ps_phylum_rel <- transform_sample_counts(ps_phylum, function(x) x / sum(x))
phylum_df <- psmelt(ps_phylum_rel)

# Alphabetical order for legend
phylum_df <- phylum_df %>%
  arrange(River, Sample) %>%
  mutate(Phylum = factor(Phylum, levels = sort(unique(Phylum))),
         Sample = factor(Sample, levels = unique(Sample)))

# Colour palette
palette_colors <- scales::hue_pal()(length(levels(phylum_df$Phylum)))

# Plot grouped by River
p_phylum <- ggplot(phylum_df, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = palette_colors,
                    limits = levels(phylum_df$Phylum),
                    drop = FALSE) +
  facet_wrap(~ River, scales = "free_x", nrow = 1) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 16, face = "bold"),
    panel.grid.major.x = element_blank()
  ) +
  labs(x = "Sample", y = "Relative Abundance (%)",
       title = "Relative Abundance at Phylum Level") +
  guides(fill = guide_legend(ncol = 2))

# Save plot
ggsave(file.path(path, "Result/Phylum_Relative_Abundance_GroupedByRiver.png"),
       plot = p_phylum, width = 20, height = 10, dpi = 300)

# Export Phylum-level table
phylum_table <- phylum_df %>%
  group_by(Sample, Phylum) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  mutate(Abundance = round(Abundance * 100, 2)) %>%
  pivot_wider(names_from = Phylum, values_from = Abundance, values_fill = 0)

write_xlsx(phylum_table, file.path(path, "Result/Phylum_Relative_Abundance.xlsx"))

# --- Family Level ---
# Prepare Family-level data
ps_family <- tax_glom(ps_filtered, taxrank = "Family")
ps_family_rel <- transform_sample_counts(ps_family, function(x) x / sum(x))
family_df <- psmelt(ps_family_rel)

# Top 30 families
top_families <- family_df %>%
  group_by(Family) %>%
  summarise(TotalAbundance = sum(Abundance)) %>%
  arrange(desc(TotalAbundance)) %>%
  slice(1:30) %>%
  pull(Family) %>%
  as.character() %>%
  setdiff("Other")

# Recode and factor (Other last)
family_df <- family_df %>%
  mutate(Family = ifelse(Family %in% top_families, Family, "Other"),
         Family = factor(Family, levels = c(sort(top_families), "Other"))) %>%
  arrange(River, Sample) %>%
  mutate(Sample = factor(Sample, levels = unique(Sample)))

# Colour palette
palette_colors <- scales::hue_pal()(length(levels(family_df$Family)))

# Plot grouped by River
p_family <- ggplot(family_df, aes(x = Sample, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = palette_colors,
                    limits = levels(family_df$Family),
                    drop = FALSE) +
  facet_wrap(~ River, scales = "free_x", nrow = 1) +  # Group by River
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 16, face = "bold"), # River labels
    panel.grid.major.x = element_blank()
  ) +
  labs(x = "Sample", y = "Relative Abundance (%)",
       title = "Relative Abundance at Family Level (Top 30 + Other)") +
  guides(fill = guide_legend(ncol = 2))

# Save plot
ggsave(file.path(path, "Result/Family_Relative_Abundance_Top30_OtherLast_GroupedByRiver.png"),
       plot = p_family, width = 20, height = 10, dpi = 300)

# Export Family-level table
family_table <- family_df %>% group_by(Sample, Family) %>% summarise(Abundance = sum(Abundance), .groups = "drop") %>% mutate(Abundance = round(Abundance * 100, 2)) %>% pivot_wider(names_from = Family, values_from = Abundance, values_fill = 0)
write_xlsx(family_table, file.path(path, "Result/Family_Relative_Abundance_Top30.xlsx"))

# --- Genus Level ---
# Prepare Genus-level data
ps_genus <- tax_glom(ps_filtered, taxrank = "Genus")
ps_genus_rel <- transform_sample_counts(ps_genus, function(x) x / sum(x))
genus_df <- psmelt(ps_genus_rel)

# Top 30 genera
top_genera <- genus_df %>%
  group_by(Genus) %>%
  summarise(TotalAbundance = sum(Abundance)) %>%
  arrange(desc(TotalAbundance)) %>%
  slice(1:30) %>%
  pull(Genus) %>%
  as.character() %>%
  setdiff("Other")

# Recode and factor (Other last)
genus_df <- genus_df %>%
  mutate(Genus = ifelse(Genus %in% top_genera, Genus, "Other"),
         Genus = factor(Genus, levels = c(sort(top_genera), "Other"))) %>%
  arrange(River, Sample) %>%
  mutate(Sample = factor(Sample, levels = unique(Sample)))

# Colour palette
palette_colors <- scales::hue_pal()(length(levels(genus_df$Genus)))

# Plot grouped by River
p_genus <- ggplot(genus_df, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = palette_colors,
                    limits = levels(genus_df$Genus),
                    drop = FALSE) +
  facet_wrap(~ River, scales = "free_x", nrow = 1) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 16, face = "bold"),
    panel.grid.major.x = element_blank()
  ) +
  labs(x = "Sample", y = "Relative Abundance (%)",
       title = "Relative Abundance at Genus Level (Top 30 + Other)") +
  guides(fill = guide_legend(ncol = 2))

# Save plot
ggsave(file.path(path, "Result/Genus_Relative_Abundance_Top30_OtherLast_GroupedByRiver.png"),
       plot = p_genus, width = 20, height = 10, dpi = 300)

#Export Genus-level table
genus_table <- genus_df %>% group_by(Sample, Genus) %>% summarise(Abundance = sum(Abundance), .groups = "drop") %>% mutate(Abundance = round(Abundance * 100, 2)) %>% pivot_wider(names_from = Genus, values_from = Abundance, values_fill = 0)
write_xlsx(genus_table, file.path(path, "Result/Genus_Relative_Abundance_Top30.xlsx"))
