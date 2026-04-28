Microbiome Analysis (16S rRNA)

This repository contains microbiome data analysis using river water samples. The analysis follows a standard 16S rRNA amplicon sequencing workflow in R, from sequence processing to community-level analysis.

Objectives:
1. To process raw 16S rRNA sequencing data into ASVs
2. To characterize microbial community composition
3. To analyse alpha and beta diversity across river water 

Workflow:
1. Sequence Processing (DADA2):
    a. Quality filtering and trimming
    b. Error model learning
    c. Dereplication and ASV inference
    d. Chimera removal
3. Taxonomic Assignment:
    a. Classification using the SILVA reference database
4. Phyloseq Analysis:
    a. Integration of ASV table, taxonomy, and metadata
    b. Data transformation to relative abundance
    c. Taxonomic aggregation (e.g., phylum/genus level)

Analyses Performed:
    1. Relative abundance bar plots
    2. Alpha diversity analysis
    3. Beta diversity analysis (ordination: PCoA/NMDS)
    4. Comparison of microbial composition across rivers

Requirements:
    1. R (≥ 4.x)
    2. Packages:
        a. dada2
        b. phyloseq
        c. ggplot2
        d. dplyr
        e. scales
        f. tidyr
        g. openxlsx
        h. writexl
        i. vegan
        j. pheatmap
        l. multcomp

Notes:
    1. This analysis was conducted as part of a microbiome data analysis workshop
    2. The dataset is used for training and demonstration purposes

Author:
    Nur Azzah Binti Osman
