# Multi-trait Finlay-Wilkinson (MT-FW) regression

This repository contains the R code and data accompanying the manuscript:

> **Dissecting the genetic architecture of multi-trait plasticity across crops.**

The pipeline estimates genotype **plasticity slopes** with single-trait and
multi-trait Finlay–Wilkinson (FW) regression — the latter built on non-negative matrix
decompositions, namely the classical non-negative matrix factororization (NMF) and the non-negative 
principal component analysis (PCA) of the multi-trait phenotype matrix;
and then determines the genetic architecture of those slopes in GWAS using FarmCPU from
GAPIT 3.

---

## Table of contents

- [Repository structure](#repository-structure)
- [Data files](#data-files)
- [Software requirements](#software-requirements)
- [How to reproduce the analysis](#how-to-reproduce-the-analysis)
- [Script descriptions](#script-descriptions)
- [License](#license)

---

## Repository structure

```
.
├── README.md
├── LICENSE
├── .gitignore
├── data/
│   ├── Maize_Pheno_Data_4Env.csv          # Maize: agronomic traits + NIRS, 4 environments
│   ├── Coffee_Pheno_Data_3Env.csv         # Coffee: agronomic + photosynthetic, 3 environments
│   ├── my_combined_betas_Maize.csv        # Combined ST + MT FW slopes (GWAS phenotype input)
│   ├── SNP60000_Common146_Maize.hmp.txt   # Maize GBS SNPs (HapMap format, 146 genotypes)
│   └── table_GWAS_Results_snpbiall.csv    # Summary table of all GWAS associations
│
└── R/
    ├── 01_finlay_wilkinson/
    │   └── FW_ST_MT.R                      # ST-FW, NMF / nnPCA decomposition, MT-FW
    │
    └── 02_gwas/
        ├── GWAS.R                          # FarmCPU (GAPIT 3); effective number of tests
        └── Modified_table_GWAS_Results.R  # 3-threshold significance table function
```

> **Note on the GAPIT functions file.** `GWAS.R` sources the GAPIT 3 functions at
> run time from a local `gapit_functions.txt` (a bug-fixed version of the upstream
> file). Because this file is distributed by the GAPIT authors, download the
> current version from <https://zzlab.net/GAPIT/> and place it next to `GWAS.R`
> rather than redistributing it here.

---

## Data files

| File | Description |
|------|-------------|
| `Maize_Pheno_Data_4Env.csv` | Maize hybrid panel: agronomic traits + selected NIRS spectral bands, 4 environments (146 genotypes) |
| `Coffee_Pheno_Data_3Env.csv` | Three-way hybrid coffee population: agronomic and photosynthetic traits, 3 environments |
| `my_combined_betas_Maize.csv` | Combined single-trait and multi-trait FW slopes; first column = genotype; used directly as the GWAS phenotype input |
| `SNP60000_Common146_Maize.hmp.txt` | Maize GBS SNP genotypes in HapMap format (read with `header = FALSE`) |
| `table_GWAS_Results_snpbiall.csv` | Master summary of all GWAS associations |

---

## Software requirements

- **R** ≥ 4.2.0

CRAN packages used by the two scripts:

```r
install.packages(c(
  # Data wrangling / plotting
  "tidyverse", "dplyr", "tidyr", "purrr", "readr",
  "ggplot2", "pheatmap", "RColorBrewer", "reshape2",
  "gridExtra", "gplots", "viridis", "conflicted", "Hmisc",
  # Finlay–Wilkinson and non-negative decomposition
  "FW", "nsprcomp", "RcppML",
  # GWAS post-processing
  "data.table", "poolr"
))
```

GAPIT 3 is sourced at run time (not installed from CRAN):

```r
source("https://zzlab.net/GAPIT/GAPIT.library.R")
source("gapit_functions.txt")   # local bug-fixed copy, see note above
```

`gwaspr` (GWAS post-processing helpers used in `GWAS.R`) is on GitHub:

```r
# install.packages("devtools")
devtools::install_github("derekmichaelwright/gwaspr")
```

---

## How to reproduce the analysis

The scripts use relative paths; set the working directory to each stage folder
(or adjust paths to the repository root).

### Stage 1 — Finlay–Wilkinson slopes and meta-traits

```r
setwd("R/01_finlay_wilkinson")
source("FW_ST_MT.R")
```

Produces single-trait slopes, the non-negative decompositions (NMF, nnPCA) with
their score and loading matrices, and the multi-trait slopes.

### Stage 2 — GWAS

```r
setwd("R/02_gwas")
# Ensure gapit_functions.txt is present in this folder
source("GWAS.R")
```

Runs FarmCPU for every slope, computes the effective number of independent tests
(M_eff), and orders the result files. The
`Modified_table_GWAS_Results.R` function builds the significance table with three
thresholds (Bonferroni, effective-tests, exploratory).

---

## Script descriptions

### `R/01_finlay_wilkinson/FW_ST_MT.R`

Finlay–Wilkinson workflow:
1. Environment-wise min–max scaling of the phenotype matrix.
2. **Single-trait FW** (`FW::FW`, OLS)  one slope per genotype per trait.
3. **Non-negative decomposition** of the scaled multi-trait matrix: non-negative
   PCA (`nsprcomp`) and NMF (`RcppML`); rank chosen by cumulative variance explained.
4. Reshaping of meta-trait scores into a genotype × meta-trait × environment array.
5. **Multi-trait FW** on each meta-trait → meta-trait plasticity slopes.
6. Correlation and trait–meta-trait concordance analyses (single-trait vs
   meta-trait slopes).

### `R/02_gwas/GWAS.R`

GWAS on the plasticity slopes:
- Reads the combined slope matrix (`my_combined_betas_Maize.csv`) as phenotype and
  the maize GBS HapMap (`SNP60000_Common146_Maize.hmp.txt`) as genotype.
- Derives PCA and kinship once, for consistency across all slopes.
- Runs **FarmCPU** (GAPIT 3) for every slope.
- Computes the **effective number of tests**  for the effective-tests threshold.

### `R/02_gwas/Modified_table_GWAS_Results.R`

A modified GAPIT/`gwaspr` table function that classifies each SNP under a
**three-threshold** scheme (Significant / Suggestive / Exploratory), enabling the
exploratory threshold (e.g. −log₁₀P = 4.0) used in the paper.

---

## License

Released under the MIT License — see [`LICENSE`](LICENSE).

The phenotype and genotype data are from previous publications and are provided here for reproducibility only. 
The GAPIT functions file is the property of the GAPIT
authors and is **not** redistributed here; download it from
<https://zzlab.net/GAPIT/>.

---

