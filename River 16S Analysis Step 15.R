River_16S_Analysis.R (with Unfiltered + Filtered branches)
Step 15
###############################################################################
# River Microbiome Analysis (16S rRNA)
# Author: Azzah
# Updated: 2025-11-13
# Description: Full DADA2 + Phyloseq pipeline with unfiltered alpha diversity
#              and filtered beta diversity / relative abundance / heatmap
###########################################################################
####
# 15. ALPHA DIVERSITY (Unfiltered branch – preserves rare taxa)
###############################################################################
set.seed(123)
rare_depth <- min(sample_sums(ps_unfiltered))
ps_rarefied <- rarefy_even_depth(ps_unfiltered, sample.size = rare_depth, rngseed = 123, verbose = FALSE)

p_alpha <- plot_richness(ps_rarefied, x = "River", measures = c("Observed", "Shannon", "Simpson")) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
  labs(title = "Alpha Diversity (Unfiltered - Rarefied)")

# Ensure Result directory exists before saving
if (!dir.exists(file.path(path, "Result"))) dir.create(file.path(path, "Result"), recursive = TRUE)

# Save rarefied alpha diversity plot
ggsave(file.path(path, "Result/Alpha_Diversity_Rarefied.png"),
       plot = p_alpha, width = 10, height = 8, dpi = 300)

# Save unrarefied alpha diversity metrics
alpha_df <- estimate_richness(ps_unfiltered)
write.csv(alpha_df, file.path(path, "Result/Alpha_Diversity_Unrarefied.csv"),
          row.names = TRUE)

# ANOVA on Alpha Diversity Metrics
#if not already installed
install.packages("--")
# Load required libraries
library(stats)
library(multcomp)

# Read alpha diversity table
alpha_df <- read.csv(file.path("C:/Users/erino/Desktop/13-11-2025 River R/Result/Alpha_Diversity_Unrarefied.csv"),
                     header = TRUE, row.names = 1)

# Add River column from row names (assuming SampleID format includes River)
alpha_df$SampleID <- rownames(alpha_df)
alpha_df$River <- sub("_.*", "", alpha_df$SampleID)

# Perform ANOVA for each metric
anova_results <- list()
tukey_results <- list()

metrics <- c("Observed", "Shannon", "Simpson")
for (metric in metrics) {
  formula <- as.formula(paste(metric, "~ River"))
  aov_model <- aov(formula, data = alpha_df)
  anova_results[[metric]] <- summary(aov_model)

  # Tukey HSD if significant
  tukey <- TukeyHSD(aov_model)
  tukey_results[[metric]] <- tukey
}

# Save ANOVA and Tukey results to CSV
anova_out <- file.path("C:/Users/erino/Desktop/13-11-2025 River R/Result", "ANOVA_AlphaDiversity_Results.csv")
tukey_out <- file.path("C:/Users/erino/Desktop/13-11-2025 River R/Result", "TukeyHSD_AlphaDiversity_Results.csv")

# Write ANOVA results
sink(anova_out)
cat("ANOVA Results for Alpha Diversity Metrics

")
for (metric in metrics) {
  cat("Metric:", metric, "
")
  print(anova_results[[metric]])
  cat("
")
}
sink()

# Write Tukey results
sink(tukey_out)
cat("Tukey HSD Post-hoc Results

")
for (metric in metrics) {
  cat("Metric:", metric, "
")
  print(tukey_results[[metric]])
  cat("
")
}
sink()

cat("ANOVA and Tukey HSD results saved to:
", anova_out, "
", tukey_out, "
")
