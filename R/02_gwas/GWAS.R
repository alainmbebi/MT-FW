library(nlme)
library(sjstats)
library(lme4)
library(RColorBrewer)
#library(dplyr)
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
library(BGLR)
library(lessR)
library(devtools)
library(webr)
library(tidyverse) 
library(pheatmap) 
library(RColorBrewer) 
#library(tidyr)
library(purrr)
library(readr)
library(FW) #Package to perform FWR

# =========================================================
# =========================================================
# =========================================================
#Source function and gapit library
#### HapMap data ####

source("http:/zzlab.net/GAPIT/GAPIT.library.R") 
#source("http:/zzlab.net/GAPIT/gapit_functions.txt")
#try new function below as the function above failed in gwas. 
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
  Random.model = T,  # Optional: use if GAPIT returns an error
  Phenotype.View = T ,# Optional: use if GAPIT returns an error
  file.output = TRUE
)

### ---------------------------------------------------------
#Plot the derived data to determine the optimal PCs
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

#Compute PCA (population structure)
PCA_matrix <- myGAPIT_PCA$PCA #we will reuse this PCA (important for consistency across traits)
#Compute Kinship (K matrix)
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
threshold <- alpha / Meff_total ####0.00002343018
eff.threshold = -log10(threshold) ###4.630224
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

#FIGURES FOR THE ARTICLE GWAS PART
# ============================================================
# ============================================================
# LD DECAY ANALYSIS - Zea mays GBS (TASSEL single-file input)
# Remington 2001 HW-adjusted model + LOESS + rolling mean
# Per-chromosome and genome-wide average
# ============================================================

rm(list = ls())
gc()

# Install and load packages
pkgs_cran <- c("data.table", "ggplot2", "minpack.lm", "zoo",
               "patchwork", "viridis", "scales", "dplyr")

for (pkg in pkgs_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

library(data.table)
library(ggplot2)
library(minpack.lm)
library(zoo)
library(patchwork)
library(viridis)
library(scales)
library(dplyr)

# ============================================================
# SETTINGS
# ============================================================

input_file  <- "LD_Maize_GBS.txt" #DErived from TASSEL 5
output_dir  <- "LD_results"
sample_prop <- 0.1        # proportion to sample for LOESS/NLS
uncorr_dist <- 5e5        # 500 kb threshold for background LD

dir.create(output_dir, showWarnings = FALSE)

# ============================================================
# FUNCTIONS
# ============================================================

LD_decay_HW_adj <- function(d, c, n) {
  C <- c * d
  ((10 + C) / ((2 + C) * (11 + C))) *
    (1 + ((3 + C) * (12 + 12 * C + C^2)) /
       (n * (2 + C) * (11 + C)))
}

LD_cut_HW_adj <- function(d, c, r, n) {
  LD_decay_HW_adj(d, c, n) - r
}

# ============================================================
# READ AND PREPARE DATA
# ============================================================

cat("Reading LD file (large file, please wait)...\n")

LD_all <- fread(input_file, header = TRUE, sep = "\t")
colnames(LD_all) <- gsub("\\^", "", colnames(LD_all))

# Numeric conversion on key columns
LD_all[, Dist_bp := as.numeric(Dist_bp)]
LD_all[, R2      := as.numeric(R2)]
LD_all[, pDiseq  := as.numeric(pDiseq)]
LD_all[, Locus1  := as.integer(Locus1)]

# Verify conversion
cat("Column classes after conversion:\n")
cat(sprintf("  Dist_bp : %s\n", class(LD_all$Dist_bp)))
cat(sprintf("  R2      : %s\n", class(LD_all$R2)))
cat(sprintf("  pDiseq  : %s\n", class(LD_all$pDiseq)))

cat(sprintf("Loaded %d SNP pairs\n", nrow(LD_all)))
cat("Columns:", paste(names(LD_all), collapse = ", "), "\n")

stopifnot("R2"      %in% colnames(LD_all))
stopifnot("Dist_bp" %in% colnames(LD_all))
stopifnot("Locus1"  %in% colnames(LD_all))
stopifnot("pDiseq"  %in% colnames(LD_all))

LD_all <- LD_all[
  !is.na(R2) & !is.na(Dist_bp) &
    R2 >= 0 & R2 <= 1 &
    Dist_bp > 0 &
    pDiseq < 0.05 &
    Locus1 %in% 1:10
]

LD_all[, Chr := as.integer(Locus1)]
setorder(LD_all, Chr, Dist_bp)

cat(sprintf("After filtering: %d SNP pairs\n", nrow(LD_all)))
cat("Pairs per chromosome:\n")
print(table(LD_all$Chr))

# ============================================================
# COLOR PALETTE AND THEME
# ============================================================

chr_colors <- viridis(10, option = "turbo")
names(chr_colors) <- paste0("Chr", 1:10)

pub_theme <- theme_classic(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 11, hjust = 0),
    plot.subtitle    = element_text(size = 8, color = "grey40", hjust = 0),
    axis.text        = element_text(size = 9, color = "black"),
    axis.title       = element_text(size = 10, face = "bold"),
    panel.grid.major = element_line(color = "grey93", linewidth = 0.3),
    legend.position  = "none",
    plot.margin      = margin(6, 8, 6, 6)
  )

# ============================================================
# PER-CHROMOSOME ANALYSIS
# ============================================================

summary_list <- list()
plot_list    <- list()

outfile <- file.path(output_dir, "LD_summary.txt")
cat("Chr\tBG_LD\tSampleSize\tLD_decay_bp\tLD_decay_kb\tMethod\n",
    file = outfile)

for (chr in 1:10) {
  
  cat(sprintf("\n--- Processing Chromosome %d ---\n", chr))
  
  data_chr  <- LD_all[Chr == chr, .(Dist = Dist_bp, R2 = R2)]
  data_chr  <- na.omit(data_chr)
  setorder(data_chr, Dist)
  data_full <- as.data.frame(data_chr)
  
  if (nrow(data_full) < 100) {
    cat(sprintf("  Chr %d: too few pairs (%d), skipping.\n", chr, nrow(data_full)))
    next
  }
  
  # DIAGNOSTICS CHECK
  data_uncorr <- data_full[data_full$Dist > uncorr_dist, ]
  cat(sprintf("  Total pairs        : %d\n", nrow(data_full)))
  cat(sprintf("  Pairs > %.0f kb    : %d\n", uncorr_dist / 1000, nrow(data_uncorr)))
  cat(sprintf("  Max distance       : %.0f bp (%.2f kb)\n",
              max(data_full$Dist), max(data_full$Dist) / 1000))
  
  #  BACKGROUND LD 
  background_LD <- if (nrow(data_uncorr) >= 10) {
    quantile(data_uncorr$R2, 0.95, na.rm = TRUE)
  } else {
    quantile(
      data_full$R2[data_full$Dist > quantile(data_full$Dist, 0.8)],
      0.95, na.rm = TRUE
    )
  }
  cat(sprintf("  Background LD (95th pct): %.4f\n", background_LD))
  
  #  SAMPLE 
  set.seed(123)
  sam_size <- min(max(floor(nrow(data_full) * sample_prop), 500), nrow(data_full))
  data_sam <- na.omit(data_full[sample(nrow(data_full), sam_size), ])
  
  #  LOESS FIT 
  fit_loess  <- loess(R2 ~ Dist, data = data_sam, degree = 2, span = 0.1)
  loess_x    <- sort(data_sam$Dist)
  loess_y    <- predict(fit_loess, newdata = data.frame(Dist = loess_x))
  data_loess <- data.frame(Dist = loess_x, R2 = pmax(0, loess_y))
  data_loess <- data_loess[order(data_loess$Dist), ]
  
  #  ROLLING MEAN 
  roll_mean  <- zoo::rollapply(data_full$R2, width = 100,
                               FUN = mean, fill = NA, align = "center")
  rolling_df <- data.frame(Dist = data_full$Dist, R2 = roll_mean)
  
  #  HW-ADJUSTED MODEL (Remington 2001) 
  fit_HW_adj <- try(
    nlsLM(
      R2 ~ LD_decay_HW_adj(Dist, c, sam_size),
      data    = data_sam,
      start   = list(c = 0.1),
      lower   = 0,
      control = nls.lm.control(maxiter = 200)
    ),
    silent = TRUE
  )
  
  c_val <- NA; LD_hw <- NA; hw_df <- NULL
  
  if (!inherits(fit_HW_adj, "try-error")) {
    c_val <- coef(fit_HW_adj)["c"]
    d_seq <- seq(1, max(data_full$Dist), length.out = 1000)
    hw_df <- data.frame(Dist = d_seq,
                        R2   = LD_decay_HW_adj(d_seq, c_val, sam_size))
    root  <- try(
      uniroot(LD_cut_HW_adj,
              interval = c(1, max(data_full$Dist)),
              c = c_val, r = background_LD, n = sam_size)$root,
      silent = TRUE
    )
    if (!inherits(root, "try-error")) LD_hw <- root
  } else {
    cat("  HW_adj model failed - using LOESS only\n")
  }
  
  # --- LOESS DECAY DISTANCE ---
  valid_idx <- which(data_loess$R2 > background_LD)
  LD_loess  <- if (length(valid_idx) > 0) {
    max(data_loess$Dist[valid_idx], na.rm = TRUE)
  } else NA
  
  # --- FINAL DECAY DISTANCE (HW_adj > LOESS > approx > fallback) ---
  if (!is.na(LD_hw) && is.numeric(LD_hw)) {
    LD_final    <- LD_hw
    method_used <- "HW_adj"
  } else if (!is.na(LD_loess) && is.numeric(LD_loess)) {
    LD_final    <- LD_loess
    method_used <- "LOESS"
  } else {
    idx <- which.min(abs(data_loess$R2 - background_LD))
    if (length(idx) > 0 && !is.na(data_loess$Dist[idx])) {
      LD_final    <- data_loess$Dist[idx]
      method_used <- "approx"
    } else {
      LD_final    <- max(data_full$Dist)
      method_used <- "fallback_max"
      cat("  WARNING: Could not determine LD decay - using max distance\n")
    }
  }
  
  LD_final    <- as.numeric(LD_final)
  LD_final_kb <- round(LD_final / 1000, 2)
  cat(sprintf("  LD decay: %.0f bp (%.2f kb) [%s]\n",
              LD_final, LD_final_kb, method_used))
  
  #  SAVE 
  cat(sprintf("%d\t%.4f\t%d\t%.0f\t%.2f\t%s\n",
              chr, background_LD, sam_size,
              LD_final, LD_final_kb, method_used),
      file = outfile, append = TRUE)
  
  write.csv(
    data.frame(Chr = chr, BG_LD = background_LD,
               LD_decay_bp = LD_final, LD_decay_kb = LD_final_kb,
               Method = method_used),
    file      = file.path(output_dir, paste0("LD_decay_chr", chr, ".csv")),
    row.names = FALSE
  )
  
  summary_list[[chr]] <- data.frame(
    Chr         = paste0("Chr", chr),
    BG_LD       = round(background_LD, 4),
    SampleSize  = sam_size,
    LD_decay_bp = round(LD_final, 0),
    LD_decay_kb = LD_final_kb,
    Method      = method_used
  )
  
  #  PER-CHROMOSOME PLOT 
  chr_label <- paste0("Chr", chr)
  col_chr   <- chr_colors[chr_label]
  
  set.seed(42)
  plot_pts <- data_full[sample(nrow(data_full), min(5000, nrow(data_full))), ]
  
  p <- ggplot() +
    geom_point(data = plot_pts,
               aes(x = Dist / 1000, y = R2),
               color = col_chr, size = 0.25, alpha = 0.25, shape = 16) +
    geom_line(data = rolling_df,
              aes(x = Dist / 1000, y = R2),
              color = "steelblue", linewidth = 0.5, na.rm = TRUE) +
    geom_line(data = data_loess,
              aes(x = Dist / 1000, y = R2),
              color = "red", linewidth = 0.8) +
    geom_hline(yintercept = background_LD, color = "purple",
               linetype = "dashed", linewidth = 0.5) +
    geom_vline(xintercept = LD_final / 1000, color = "black",
               linetype = "dashed", linewidth = 0.5) +
    annotate("text",
             x = LD_final / 1000, y = 0.92,
             label    = paste0(LD_final_kb, " kb"),
             hjust    = -0.1, size = 2.8,
             color    = "black", fontface = "bold") +
    scale_y_continuous(limits = c(0, 1),
                       breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_x_continuous(labels = comma_format()) +
    labs(
      title = paste0("Chr ", chr,
                     "  |  decay: ", LD_final_kb, " kb",
                     "  |  bg R2: ", round(background_LD, 3)),
                    # "  [", method_used, "]"), #uncomment here to have this on the title
      x = "Distance (kb)",
      y = expression(R^2)
    ) +
    pub_theme
  
  if (!is.null(hw_df)) {
    p <- p + geom_line(data = hw_df, aes(x = Dist / 1000, y = R2),
                       color = "darkgreen", linetype = "dashed", linewidth = 0.7)
  }
  
  plot_list[[chr]] <- p
}

# ============================================================
# SUMMARY TABLE
# ============================================================

summary_df <- do.call(rbind, summary_list)
rownames(summary_df) <- NULL

cat("\n========== LD DECAY SUMMARY ==========\n")
print(summary_df)

avg_decay_kb <- round(mean(summary_df$LD_decay_kb, na.rm = TRUE), 2)
avg_decay_bp <- round(mean(summary_df$LD_decay_bp, na.rm = TRUE), 0)
cat(sprintf("\nGenome-wide average LD decay: %.0f bp (%.2f kb)\n",
            avg_decay_bp, avg_decay_kb))

write.csv(summary_df,
          file.path(output_dir, "LD_decay_all_chromosomes.csv"),
          row.names = FALSE)

# ============================================================
# FIGURE 1: PER-CHROMOSOME PANEL (2 rows x 5 cols)
# ============================================================

cat("\nGenerating per-chromosome figure...\n")

fig_chr <- wrap_plots(plot_list, ncol = 5, nrow = 2) +
  plot_annotation(
    title    = " ",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 9, color = "grey40")
    )
  )

ggsave(file.path(output_dir, "Fig1_LD_decay_per_chromosome.pdf"),
       fig_chr, width = 20, height = 10, units = "in", dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Fig1_LD_decay_per_chromosome.png"),
       fig_chr, width = 20, height = 10, units = "in", dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Fig1_LD_decay_per_chromosome.eps"),
       fig_chr, width = 20, height = 10, units = "in", dpi = 300, bg = "white")
cat("  Saved: Fig1_LD_decay_per_chromosome.pdf/png/eps\n")

# ============================================================
# FIGURE 2: GENOME-WIDE AVERAGE LD DECAY
# ============================================================

cat("Generating genome-wide average figure...\n")

data_avg        <- as.data.frame(LD_all[, .(Dist = Dist_bp, R2, Chr)])
data_avg_sorted <- data_avg[order(data_avg$Dist), ]

set.seed(42)
sam_avg      <- min(max(floor(nrow(data_avg) * sample_prop), 1000), nrow(data_avg))
data_sam_avg <- data_avg[sample(nrow(data_avg), sam_avg), ]

roll_avg       <- zoo::rollapply(data_avg_sorted$R2, width = 500,
                                 FUN = mean, fill = NA, align = "center")
rolling_avg_df <- data.frame(Dist = data_avg_sorted$Dist, R2 = roll_avg)

fit_loess_avg  <- loess(R2 ~ Dist, data = data_sam_avg, degree = 2, span = 0.1)
loess_avg_x    <- sort(data_sam_avg$Dist)
loess_avg_y    <- predict(fit_loess_avg, newdata = data.frame(Dist = loess_avg_x))
data_loess_avg <- data.frame(Dist = loess_avg_x, R2 = pmax(0, loess_avg_y))
data_loess_avg <- data_loess_avg[order(data_loess_avg$Dist), ]

data_uncorr_avg <- data_avg[data_avg$Dist > uncorr_dist, ]
bg_avg <- if (nrow(data_uncorr_avg) >= 10) {
  quantile(data_uncorr_avg$R2, 0.95, na.rm = TRUE)
} else {
  quantile(data_avg$R2[data_avg$Dist > quantile(data_avg$Dist, 0.8)],
           0.95, na.rm = TRUE)
}
cat(sprintf("  Genome-wide background LD: %.4f\n", bg_avg))

fit_HW_avg <- try(
  nlsLM(R2 ~ LD_decay_HW_adj(Dist, c, sam_avg),
        data = data_sam_avg, start = list(c = 0.1), lower = 0,
        control = nls.lm.control(maxiter = 200)),
  silent = TRUE
)

hw_avg_df <- NULL; LD_hw_avg <- NA; c_val_avg <- NA

if (!inherits(fit_HW_avg, "try-error")) {
  c_val_avg <- coef(fit_HW_avg)["c"]
  d_seq     <- seq(1, max(data_avg$Dist), length.out = 1000)
  hw_avg_df <- data.frame(Dist = d_seq,
                          R2   = LD_decay_HW_adj(d_seq, c_val_avg, sam_avg))
  root_avg  <- try(
    uniroot(LD_cut_HW_adj,
            interval = c(1, max(data_avg$Dist)),
            c = c_val_avg, r = bg_avg, n = sam_avg)$root,
    silent = TRUE
  )
  if (!inherits(root_avg, "try-error")) LD_hw_avg <- root_avg
}

valid_avg    <- which(data_loess_avg$R2 > bg_avg)
LD_loess_avg <- if (length(valid_avg) > 0) {
  max(data_loess_avg$Dist[valid_avg], na.rm = TRUE)
} else NA

if (!is.na(LD_hw_avg) && is.numeric(LD_hw_avg)) {
  LD_avg_final <- LD_hw_avg; method_avg <- "HW_adj"
} else if (!is.na(LD_loess_avg) && is.numeric(LD_loess_avg)) {
  LD_avg_final <- LD_loess_avg; method_avg <- "LOESS"
} else {
  LD_avg_final <- avg_decay_bp; method_avg <- "mean_of_chr"
}

LD_avg_kb <- round(as.numeric(LD_avg_final) / 1000, 2)
cat(sprintf("  Genome-wide avg LD decay: %.0f bp (%.2f kb) [%s]\n",
            LD_avg_final, LD_avg_kb, method_avg))

set.seed(99)
plot_pts_avg          <- data_avg[sample(nrow(data_avg), min(20000, nrow(data_avg))), ]
plot_pts_avg$ChrLabel <- paste0("Chr", plot_pts_avg$Chr)

fig_avg <- ggplot() +
  geom_point(data = plot_pts_avg,
             aes(x = Dist / 1000, y = R2, color = ChrLabel),
             size = 0.2, alpha = 0.15, shape = 16) +
  geom_line(data = rolling_avg_df,
            aes(x = Dist / 1000, y = R2),
            color = "steelblue", linewidth = 0.5, na.rm = TRUE) +
  geom_line(data = data_loess_avg,
            aes(x = Dist / 1000, y = R2),
            color = "red", linewidth = 1.0) +
  geom_hline(yintercept = bg_avg, color = "purple",
             linetype = "dashed", linewidth = 0.6) +
  geom_vline(xintercept = LD_avg_final / 1000, color = "black",
             linetype = "dashed", linewidth = 0.7) +
  annotate("text",
           x = LD_avg_final / 1000, y = 0.93,
           label    = paste0("Avg decay:\n", LD_avg_kb, " kb"),
           hjust    = -0.08, size = 3.8,
           color    = "black", fontface = "bold") +
  annotate("text",
           x     = max(plot_pts_avg$Dist / 1000) * 0.97,
           y     = bg_avg,
           label = paste0("Background R2 = ", round(bg_avg, 3)),
           hjust = 1, vjust = -0.5, size = 3.2,
           color = "purple", fontface = "italic") +
  scale_color_manual(values = chr_colors, name = "Chromosome") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(labels = comma_format()) +
  guides(color = guide_legend(override.aes = list(alpha = 0.9, size = 2), ncol = 2)) +
  labs(
    title    = " ",
    x = "Distance (kb)", y = expression(R^2)
  ) +
  pub_theme +
  theme(legend.position = "right")

if (!is.null(hw_avg_df)) {
  fig_avg <- fig_avg +
    geom_line(data = hw_avg_df, aes(x = Dist / 1000, y = R2),
              color = "darkgreen", linetype = "dashed", linewidth = 0.8)
}

ggsave(file.path(output_dir, "Fig2_LD_decay_average.pdf"),
       fig_avg, width = 12, height = 7, units = "in", dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Fig2_LD_decay_average.png"),
       fig_avg, width = 12, height = 7, units = "in", dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Fig2_LD_decay_average.eps"),
       fig_avg, width = 12, height = 7, units = "in", dpi = 300, bg = "white")
cat("  Saved: Fig2_LD_decay_average.pdf/png/eps\n")

# ============================================================
# FIGURE 3: SUMMARY BARPLOT
# ============================================================

cat("Generating summary barplot...\n")

summary_df$Chr <- factor(summary_df$Chr, levels = paste0("Chr", 1:10))

fig_bar <- ggplot(summary_df, aes(x = Chr, y = LD_decay_kb, fill = Chr)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.3) +
  geom_hline(yintercept = LD_avg_kb, color = "red",
             linetype = "dashed", linewidth = 0.8) +
  geom_text(aes(label = paste0(LD_decay_kb, " kb")),
            vjust = -0.4, size = 3.2, fontface = "bold", color = "grey20") +
  annotate("text",
           x = 10.4, y = LD_avg_kb,
           label    = paste0("Avg: ", LD_avg_kb, " kb"),
           hjust    = 1, vjust = -0.5,
           color    = "red", size = 3.5, fontface = "bold") +
  scale_fill_manual(values = chr_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title    = " ", 
    x = "Chromosome", y = "LD decay distance (kb)"
  ) +
  pub_theme +
  theme(legend.position = "none",
        axis.text.x     = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "Fig3_LD_decay_summary.pdf"),
       fig_bar, width = 9, height = 6, units = "in", dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Fig3_LD_decay_summary.png"),
       fig_bar, width = 9, height = 6, units = "in", dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Fig3_LD_decay_summary.eps"),
       fig_bar, width = 9, height = 6, units = "in", dpi = 300, bg = "white")
cat("  Saved: Fig3_LD_decay_summary.pdf/png/eps\n")


