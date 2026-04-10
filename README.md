# Equine_skin_scSeq
A single-cell RNA-seq dataset characterizing cellular diversity in healthy equine skin

# Overview
The skin serves as the primary barrier tissue in horses and is frequently affected by immune-mediated dermatological conditions, most notably insect bite hypersensitivity (IBH). This repository provides the first single-cell transcriptomic reference of normal equine skin, comprising 85,574 high-quality cells from two horses profiled in biological duplicate. Unsupervised clustering resolved 22 distinct cell populations spanning keratinocyte subpopulations, adnexal epithelial lineages, stromal and vascular compartments, and resident immune cells. This dataset provides a cellular and molecular framework for investigating equine skin diseases, wound repair, immune responses, and comparative skin biology.

# Data
Skin biopsies were collected from the mane region of two clinically healthy horses and processed using both automated and manual dissociation protocols. Single-cell libraries were prepared using the Singleron GEXSCOPE platform and sequenced on Illumina NovaSeq X (paired-end 150 bp). Raw data were aligned to the EquCab3.0 reference genome using the CeleScope pipeline.
Data Availability:
1. Raw FASTQ files: NCBI BioProject PRJNA1441223
2. Processed data: Figshare https://doi.org/10.6084/m9.figshare.31980441.


# Code Structure
This repository contains R scripts for:
1. Quality control and filtering
2. Doublet detection and removal
3. Cell type identification and annotation
4. Functional enrichment analysis

All analyses were performed in R version 4.5.0. Required packages are listed at the beginning of each script. Analysis directories are numbered sequentially according to the publication workflow.
# Funding
This work was supported by Formas - a Swedish Research Council for Sustainable Development, grant number 2023-01000 (S.W.). 
