library(nlme)
library(sjstats)
library(lme4)
library(RColorBrewer)
library(ggplot2)
library(conflicted)  
library(tidyverse)
conflict_prefer("collapse", "dplyr")
conflict_prefer("expand", "tidyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("lag", "dplyr")
conflict_prefer("pack", "tidyr")
conflict_prefer("unpack", "tidyr")
conflict_prefer("order", "genetics")
conflict_prefer("select", "dplyr")
library(reshape2)
library(gridExtra)
library(gplots)
library(viridis)
library(lessR)
library(devtools)
library(webr)
library(tidyverse) 
library(pheatmap) 
library(RColorBrewer) 
library(purrr)
library(readr)
library(FW) #Package to perform FWR

# =========================================================
# =========================================================
# =========================================================
#Source function and gapit library
#### HapMap data ####

source("http:/zzlab.net/GAPIT/GAPIT.library.R") 
source("gapit_functions.txt") #Updated version with bug fixed
library(gwaspr)
list_Traits(folder = "/GWAS_Results/")

### ---------------------------------------------------------
# Read in phenotype data
myY <- read.csv("my_combined_betas_Maize.csv", check.names = FALSE)
# Read in genotype data (note: header = F)
myG <- read.csv("SNP60000_Common146_Maize.hmp.txt", header = F)

### ---------------------------------------------------------
#Getting needed files PCA, Kinship...
myGAPIT_PCA <- GAPIT(
  Y = myY[,1:2],
  G = myG,
  PCA.total = 10,
  #model = c("GLM", "MLM", "MLMM", "FarmCPU", "BLINK"), 
  model = c("FarmCPU"), 
  Random.model = T, 
  Phenotype.View = T ,
  file.output = TRUE
)

### ---------------------------------------------------------
#Plot to determine the optimal PCs
# Prep data
xx <- read.csv("GAPIT.Genotype.PCA_eigenvalues.csv") %>% 
  mutate(PC = 1:146, `Eigen Values` = x,
         fill = ifelse(PC == 4, "PC", ""),
         `Percent of Variation` = 100 * `Eigen Values` / sum(`Eigen Values`) ) %>%
  slice(1:10)
# Plot
mp <- ggplot(xx,  aes(x = PC, y = `Percent of Variation`)) + 
  geom_line(linewidth = 1, alpha = 0.7) + 
  geom_point(aes(fill = fill), size = 3, pch = 21) +
  scale_fill_manual(values = c("darkgreen", "darkred")) +
  scale_x_continuous(breaks = 1:10) +
  theme_gwaspr(legend.position = "none") +
  labs(title = " ") 
# Save
ggsave("Fig_PCA_Eigen_Values.jpeg", width = 6, height = 4)

# PCA 
PCA_matrix <- myGAPIT_PCA$PCA
Kinship <- myGAPIT_PCA$KI
### ---------------------------------------------------------


# Custom GWAS this is a Q+K model
my_final_GAPIT_PCA <- GAPIT(
       Y = myY, # Phenotype data
       set.seed(1234),
       G = myG, # Genotype data
       output.numerical=TRUE,
       PCA.total = 10,
       model = c("FarmCPU"), 
       #SNP.impute = "Major",
       #SNP.MAF = 0.05,
       #CV = PCA_matrix[,1:3],
       #KI = Kinship,
       #file.output = F,
       Random.model = F,  
       Phenotype.View = T 
)

### ---------------------------------------------------------
### ---------------------------------------------------------
### ---------------------------------------------------------
#Post GWAS with normal pheno and snp biall
#Check which traits and models have ran with is.ran()
is_ran(folder = "/Gwas_Results_GAPIT/")
#Order results
order_GWAS_Results(folder = "/Gwas_Results_GAPIT/")
#Confirm if the results files are ordered by p.value with is_ordered(). 
is_ordered(folder = "/Gwas_Results_GAPIT/")

### ---------------------------------------------------------
### ---------------------------------------------------------
### ---------------------------------------------------------
### ---------------------------------------------------------
#Compute the effective number of tests M_eff
library(data.table)
library(poolr)
snp <- fread("SNP60000_Common146_Maize.hmp.txt",
             data.table = FALSE,
             check.names = FALSE)
#Check structure
str(snp[, 1:10])
dim(snp)

#Extract genotype matrix
# Remove metadata (first 4 columns)
geno <- snp[, -(1:4)]

# Convert to numeric (safe even if already numeric)
geno[] <- lapply(geno, function(x) as.numeric(as.character(x)))

# Convert to matrix
geno <- as.matrix(geno)

#Extract chromosome info
chrom <- as.numeric(snp$chrom)
#Transpose for LD computation
geno_t <- t(geno)   # individuals × SNPs
#Loop per chromosome for correlation 
chromosomes <- sort(unique(chrom))

meff_list <- list()

for(chr in chromosomes){
  
  cat("Processing chromosome:", chr, "\n")
  
  # SNP indices for this chromosome
  idx <- which(chrom == chr)
  
  geno_chr <- geno_t[, idx, drop = FALSE]
  
  # Remove monomorphic SNPs
  keep <- apply(geno_chr, 2, var, na.rm = TRUE) > 0
  geno_chr <- geno_chr[, keep, drop = FALSE]
  
  # OPTIONAL: reduce size if too large
  if(ncol(geno_chr) > 5000){
    set.seed(1)
    geno_chr <- geno_chr[, sample(ncol(geno_chr), 5000)]
  }
  
  # Compute LD matrix
  r_chr <- cor(geno_chr, use = "pairwise.complete.obs")
  
  # Compute Meff
  meff_list[[as.character(chr)]] <- meff(r_chr, method = "liji")
}

#Genome-wide Meff
Meff_total <- sum(unlist(meff_list))
Meff_total ###2134

#Significance threshold
alpha <- 0.05
threshold <- alpha / Meff_total
threshold ####0.00002343018
-log10(threshold) ###4.630224
### ---------------------------------------------------------

#Create a summary table of all significant associations
Meff_total <- 2134
alpha <- 0.05
threshold <- alpha / Meff_total 
eff.threshold = -log10(threshold) 
expl.threshold = 4.0

source("Modified_table_GWAS_Results.R")

x1 <- Modified_table_GWAS_Results(
  folder = "/Gwas_Results_GAPIT/", 
  #fnames = list_Result_Files("folder"),
  threshold      = 5.9,   # Bonferroni
  sug.threshold  = 4.6,   # Effective tests
  expl.threshold = 4.0     # Exploratory
)
# Save
write.csv(x1, "table_GWAS_Results_snpbiall.csv", row.names = F)

# =========================================================
# =========================================================
