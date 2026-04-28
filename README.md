Microbiome Analysis (16S rRNA)

This repository contains microbiome data analysis using river water samples. The analysis follows a standard 16S rRNA amplicon sequencing workflow in R, from sequence processing to community-level analysis.

Objectives
1. To process raw 16S rRNA sequencing data into ASVs
2. To characterize microbial community composition
3. To analyse alpha and beta diversity across river water 

Workflow
1. Sequence Processing (DADA2)
    Quality filtering and trimming
    Error model learning
    Dereplication and ASV inference
    Chimera removal
3. Taxonomic Assignment
    Classification using the SILVA reference database
4. Phyloseq Analysis
    Integration of ASV table, taxonomy, and metadata
    Data transformation to relative abundance
    Taxonomic aggregation (e.g., phylum/genus level)

Analyses Performed
    Relative abundance bar plots
    Alpha diversity analysis
    Beta diversity analysis (ordination: PCoA/NMDS)
    Comparison of microbial composition across rivers

Requirements
    R (≥ 4.x)
    Packages:
    dada2
    phyloseq
    ggplot2
    tidyverse
    vegan

Notes
    This analysis was conducted as part of a microbiome data analysis workshop
    The dataset is used for training and demonstration purposes

Author
    Nur Azzah Binti Osman
