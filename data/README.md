# Data

Input data for the pipeline.

| File | Description |
|------|-------------|
| `Maize_Pheno_Data_4Env.csv` | Maize hybrid panel phenotypes (agronomic + NIRS), 4 environments |
| `Coffee_Pheno_Data_3Env.csv` | Three-way hybrid coffee phenotypes, 3 environments |
| `my_combined_betas_Maize.csv` | Combined ST + MT FW slopes (GWAS phenotype input) |
| `SNP60000_Common146_Maize.hmp.txt` | Maize GBS SNPs, HapMap format (read with `header = FALSE`) |
| `table_GWAS_Results_snpbiall.csv` | Master GWAS association summary |

**Not included here:** `gapit_functions.txt` — download the current GAPIT 3
functions from <https://zzlab.net/GAPIT/> and place it next to `R/02_gwas/GWAS.R`.
