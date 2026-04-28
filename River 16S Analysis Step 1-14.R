River_16S_Analysis.R (with Unfiltered + Filtered branches)
Step 1-14
###############################################################################
# River Microbiome Analysis (16S rRNA)
# Author: Azzah
# Updated: 2025-11-13
# Description: Full DADA2 + Phyloseq pipeline with unfiltered alpha diversity
#              and filtered beta diversity / relative abundance / heatmap
###############################################################################
# -------------------------------
# 1. Set Working Directory & Load Libraries
# -------------------------------
path <- "C:/path to your file"
setwd(path)

#if not already installed
install.packages("--")

# Bioconductor installs
BiocManager::install(c("dada2","phyloseq","DECIPHER","Biostrings","ShortRead"))
# Optional helpful packages:
BiocManager::install(c("phyloseq","DESeq2")) 

#Load library
library(dada2)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(scales)
library(tidyr)
library(openxlsx)
library(writexl)
library(vegan)
library(pheatmap)
library(stats)
library(multcomp)

###############################################################################
# 2. RAW DATA INPUT
###############################################################################
# Identify forward and reverse FASTQ files
fnFs <- sort(list.files(".", pattern = "_R1_001\\.fastq\\.gz$", full.names = TRUE))
fnRs <- sort(list.files(".", pattern = "_R2_001\\.fastq\\.gz$", full.names = TRUE))

#if fail try
srcFs <- sort(list.files(path, pattern = "_R1_001.fastq.gz", full.names = TRUE))
srcRs <- sort(list.files(path, pattern = "_R2_001.fastq.gz", full.names = TRUE))
length(srcFs)
length(srcRs)

# Check pairing and extract sample names
stemF <- sub("_R1_001\\.fastq\\.gz$", "", basename(fnFs))
stemR <- sub("_R2_001\\.fastq\\.gz$", "", basename(fnRs))

stopifnot(all(stemF == stemR))
keys <- stemF
names(fnFs) <- names(fnRs) <- keys
sample_id <- sub("_S\\d+_L\\d{3}$", "", keys)

srcFs <- fnFs
srcRs <- fnRs

###############################################################################
# 3. QUALITY CHECK
###############################################################################
plotQualityProfile(srcFs[1:min(2, length(srcFs))])
plotQualityProfile(srcRs[1:min(2, length(srcRs))])

###############################################################################
# 4. FILTERING AND TRIMMING
###############################################################################
# Create filtered output directory
filtered_dir <- file.path(path, "filtered")
dir.create(filtered_dir, showWarnings = FALSE)

# Define output file paths
filtFs <- file.path(filtered_dir, paste0(sample_id, "_F_filt.fastq.gz"))
filtRs <- file.path(filtered_dir, paste0(sample_id, "_R_filt.fastq.gz"))

#if fail try
length(fnFs)
length(fnRs)
length(sample_id)
length(filtFs)
length(filtRs)
sample_id <- sapply(basename(fnFs), function(x) strsplit(x, "_")[[1]][1])
filtFs <- file.path(filtered_dir, paste0(sample_id, "_F_filt.fastq.gz"))
filtRs <- file.path(filtered_dir, paste0(sample_id, "_R_filt.fastq.gz"))
#rerun to recheck
length(fnFs)
length(fnRs)
length(filtFs)
length(filtRs)

fnFs <- srcFs
fnRs <- srcRs

# Filtering parameters
out <- filterAndTrim(fnFs, filtFs,
                    fnRs, filtRs,
                    truncLen = c(240,160),     # truncate forward at 250, reverse at 220
                    maxN = 0,                   # discard reads with ambiguous bases
                    maxEE = c(2, 2),            # expected errors: stricter for forward
                    truncQ = 2,                 # truncate reads at first quality score < 2
                    rm.phix = TRUE,             # remove phiX reads
                    compress = TRUE,
                    multithread = TRUE)

# Lenient filtering for environmental samples
out <- filterAndTrim(fnFs, filtFs,
                     fnRs, filtRs,
                     truncLen = c(0, 0),
                     maxN = 0,
                     maxEE = c(10, 10),
                     truncQ = 0,
                     rm.phix = TRUE,
                     compress = TRUE,
                     multithread = TRUE)
head(out)

###############################################################################
# 5. DEREPLICATION
###############################################################################
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)
names(derepFs) <- sample_id
names(derepRs) <- sample_id

###############################################################################
# 6. ERROR RATE LEARNING
###############################################################################
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

###############################################################################
# 7. DENOISING (ASV INFERENCE)
###############################################################################
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

###############################################################################
# 8. MERGE PAIRED READS
###############################################################################
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)

###############################################################################
# 9. SEQUENCE TABLE CONSTRUCTION
###############################################################################
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

###############################################################################
# 10. CHIMERA REMOVAL
###############################################################################
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE)
cat("Proportion of non-chimeric reads:", sum(seqtab.nochim) / sum(seqtab), "\n")

###############################################################################
# 11. TAXONOMIC ASSIGNMENT
###############################################################################
taxa <- assignTaxonomy(seqtab.nochim,
                       "C:/Users/erino/Desktop/13-11-2025 River R/silva_nr99_v138_train_set.fa.gz",
                       multithread = TRUE)

###############################################################################
# 12. EXPORT RESULTS
###############################################################################
new code
# Save both sequence tables (with unfiltered)
saveRDS(seqtab, file.path(path, "seqtab_unfiltered.rds"))      # Before chimera removal
saveRDS(seqtab.nochim, file.path(path, "seqtab_nochim.rds"))   # After chimera removal
saveRDS(taxa, file.path(path, "taxa.rds"))

# Convert to data frames
asv_df_unfiltered <- as.data.frame(seqtab)
asv_df_filtered   <- as.data.frame(seqtab.nochim)
taxa_df           <- as.data.frame(taxa)

# Create Excel workbook with both tables
output_path <- file.path(path, "dada2_results.xlsx")
wb <- createWorkbook()
addWorksheet(wb, "ASV_Unfiltered")
addWorksheet(wb, "ASV_Filtered")
addWorksheet(wb, "Taxonomy")

writeData(wb, "ASV_Unfiltered", asv_df_unfiltered)
writeData(wb, "ASV_Filtered", asv_df_filtered)
writeData(wb, "Taxonomy", taxa_df)

saveWorkbook(wb, output_path, overwrite = TRUE)
cat("Excel file saved at:", output_path, "\n")

------
old code (without unfiltered)
saveRDS(seqtab.nochim, "seqtab_nochim.rds")
saveRDS(taxa, "taxa.rds")

asv_df <- as.data.frame(seqtab.nochim)
taxa_df <- as.data.frame(taxa)
output_path <- file.path(path, "dada2_results.xlsx")

wb <- createWorkbook()
addWorksheet(wb, "ASV_Table")
addWorksheet(wb, "Taxonomy")
writeData(wb, "ASV_Table", asv_df)
writeData(wb, "Taxonomy", taxa_df)
saveWorkbook(wb, output_path, overwrite = TRUE)
cat("Excel file saved at:", output_path, "\n")

###############################################################################
# 13. PHYLOSEQ INTEGRATION
###############################################################################
metadata_raw <- read.csv(file.path(path, "metadata.csv"),
                         header = TRUE, sep = ",", stringsAsFactors = FALSE)

metadata_raw <- metadata_raw[!duplicated(metadata_raw$Sample.ID), ]
metadata_raw$Sample.ID <- gsub(" ", "_", metadata_raw$Sample.ID)
metadata_raw <- metadata_raw[metadata_raw$Sample.ID != "", ]
rownames(metadata_raw) <- metadata_raw$Sample.ID
metadata_raw$Sample.ID <- NULL
metadata <- metadata_raw

# Create unfiltered Phyloseq object (no ASV filtering, just chimera removed)
ps_unfiltered <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE),
  tax_table(taxa),
  sample_data(metadata)
)

saveRDS(ps_unfiltered, "C:/Users/erino/Desktop/13-11-2025 River R/Result/ps_unfiltered_seqtab_nochim.rds")

###############################################################################
# 14. SAMPLE FILTERING (low-read samples only)
###############################################################################
summary(sample_sums(ps_unfiltered))
hist(sample_sums(ps_unfiltered), breaks = 40, main = "Library Sizes", xlab = "Reads per sample")
min_reads <- 1000  # Adjust based on histogram
ps_unfiltered <- prune_samples(sample_sums(ps_unfiltered) >= min_reads, ps_unfiltered)
