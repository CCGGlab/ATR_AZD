This notebook contains high level information on the bioinformatics (RNA-Seq & proteomics) pipeline that was used for the analysis reported in Borenas et al. ALK signalling primes the DNA damage response sensitizing ALK-driven neuroblastoma to therapeutic ATR inhibition, PNAS 2024.  

# Environment
  
Analysis was performed in a Conda environment. See **ATR_AZD.yml** for details. **scripts/Rpacks** describes R packages that were installed independently.

# Experimental conditions

Phosphoproteomics data of CLB-BAR NB cell lines after ATR inhibition with elimusertib (BAY 1895344; 50nM) or ceralasertib (AZD; 50nM and 1µM).  Cells were synchronized and treated for 6h (DMSO control)

RNA-Seq data ALK mice after 3d ATR inhibition with 1µM ceralasertib. Data compared to previous elimusertib results (see [Szydzik et al 2021](https://www.nature.com/articles/s41467-021-27057-2)).

# Data 

- Raw mice RNA-Seq data have been deposited in Arrayexpress: E-MTAB-12961

- Processed data are available in *data/ATR_AZD_data.RData*. This file containes the following objects:
  - normalized_counts_CL: DESeq2-normalized counts from cell line RNA-Seq data (from [Szydzik et al 2021](https://www.nature.com/articles/s41467-021-27057-2))
  - normalized_counts_CL_PP: Normalized counts from cell line phosphoproteomics data
  - normalized_counts_mice: DESeq2-normalized counts from mice RNA-Seq data
  - res_diff_expr_CL: DESeq2 output from cell line RNA-Seq data (from [Szydzik et al 2021](https://www.nature.com/articles/s41467-021-27057-2))
  - res_diff_expr_CL_PP: DEP output from cell line phosphoproteomics data 
  - res_diff_expr_mice: DESeq2 output from mice RNA-Seq data
  - sample_info_CL: sample information cell line RNA-Seq data
  - sample_info_CL_PP: sample information cell line phosphoproteomics data
  - sample_info_mice: sample information cell line mice data
  - MGI_to_HGNC: MGI to HGNC mapping table
  - Genesets:
    - CGP_MGI_ls: Cellular and Genetic Perturbations genesetdownloaded from MSigDB_v2023_1_Mm
    - CP_ls: Canonical pathways geneset downloaded from MSigDB v7.4
    - CSign_MGI_ls: Cellular signatures geneset downloaded from MSigDB_v2023_1_Mm
    - DDR_ls: DNA Damage Response geneset as described by KNijnenburg et al, 2018.
    - DDR_MGI_ls: DNA Damage Response geneset (mapped to MGI) as described by KNijnenburg et al, 2018.
    - GO_ls: Gene Ontology geneset downloaded from MSigDB v7.4
    - Ha_ls: Hallmark geneset downloaded from MSigDB v7.4
    - Rea_ls: Reactome geneset downloaded from MSigDB v7.4
    - Ha_MGI_ls: Hallmark geneset downloaded from MSigDB_v2023_1_Mm
    - TFT_MGI_ls: Transcription Factor Target (GTRD) geneset downloaded from MSigDB_v2023_1_Mm
    - PDB_MGI_ls: cellular signatrue geneset downloaded from PanglaoDB

# Manuscript

The main analysis, as reported in the manuscript

## Fig. 2E: DDR comparison elimusertib and lorlatinib based on previous CL data
```{r}
source("scripts/manuscript_DDR_CL.R")
```

## Fig. 4: PP comparison elimusertib with ceralasertib

### Panels A-D: Volcanos, venn, correlation & motif enrichment
```{r}
source("scripts/manuscript_CL_PP.R")
```

### Panel E: heatmap
```{r}
source("scripts/manuscript_CL_PP_hm.R")
```

## Fig. 5G-J: RNA-Seq mice after ceralasertib treatment

### Panel G-I: volcano + GSEA
```{r}
source("scripts/manuscript_mice_RNA.R")
```

### Panel J: heatmap
```{r}
source("scripts/manuscript_mice_RNA_hm.R")
```

## Fig. 6DE: Differentiation analysis on RNA-Seq mice - old data after elimusertib treatment
```{r}
source("scripts/manuscript_mice_RNA_diff.R")
```

