# Loading relevant libraries 
library(tidyverse) 
library(pheatmap) 
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(FW) #Package to perform FWR
### ---------------------------------------------------------
### 1. Read data sets
### ---------------------------------------------------------

data_pheno_MT_FW <- read.csv(
  "Maize_Pheno_Data_4Env.csv",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

#Identify Columns
trait_cols <- colnames(data_pheno_MT_FW)[-(1:5)]

#Perform Environment-wise Min–Max Scaling
data_scaled <- data_pheno_MT_FW

data_scaled[trait_cols] <- 
  ave(
    data_pheno_MT_FW[trait_cols],
    data_pheno_MT_FW$Env,
    FUN = function(x) {
      apply(x, 2, function(col) {
        min_val <- min(col, na.rm = TRUE)
        max_val <- max(col, na.rm = TRUE)
        
        if (max_val == min_val) {
          return(rep(0, length(col)))
        } else {
          return((col - min_val) / (max_val - min_val))
        }
      })
    }
  )

#==================================================================================
### Implementation and comparative analysis
#==================================================================================

### ---------------------------------------------------------
### A. Implementation of single-trait Finlay-Wilkinson (ST-FW) regression 
### ---------------------------------------------------------

#Identify trait columns
trait_cols <- colnames(data_scaled)[-(1:5)]

STFWR_results_scaled_OLS <- list()

for (trait in trait_cols) {
  
  STFWR_results_scaled_OLS[[trait]] <- FW::FW(
    y   = data_scaled[[trait]],
    VAR = data_scaled$Pedigree,
    ENV = data_scaled$Env,
    method = "OLS"
  )
}


#==========
#Save the enrire object for later use
save(STFWR_results_scaled_OLS, file = "STFWR_results_scaled_OLS.Rdata")
#Save as matrix the slope for each genotype and traits
beta_STFW_scaled_OLS <- sapply(STFWR_results_scaled_OLS, function(res) {
  res$b
})

rownames(beta_STFW_scaled_OLS) <- unique(as.character(STFWR_results_scaled_OLS[[1]]$VAR))
write.csv(beta_STFW_scaled_OLS, "beta_STFW_scaled_OLS.csv")

#==================================================================================
### ---------------------------------------------------------
### B. Implementation of multi-trait Finlay-Wilkinson (MT-FW)regression NNPCA version
### ---------------------------------------------------------

### ---------------------------------------------------------
### ---------------------------------------------------------
#We start from the object Y_global

###----------------------
### c-Fit non-negative PCA 
###----------------------

# Identify trait columns
trait_cols <- colnames(data_scaled)[-(1:5)]

# Extract the numeric matrix

Y_global_scaled <- as.matrix(data_scaled[, trait_cols])

# Ensure non-negativity 
if (any(Y_global_scaled < 0)) stop("Negative values detected!")

library(nsprcomp)
#Important note. Two main functions are implemented in the package
#"nsprcomp" computes one principal component (PC) after
#the other. Each PA is optimized such that the corresponding PC
#has maximum additional variance not explained by the previous components
# "nscumcomp" jointly computes all PCs such that the cumulative variance is maxima
#In the current study we used the nsprcomp option

# ###----------------------
# ###----------------------
# ###----------------------
#Option non cumulative variance
set.seed(1234)
nn_pca_addi <- nsprcomp(
  Y_global_scaled,
  #ncomp = min(dim(Y_global_scaled)),   # compute full set first
  retx = TRUE,
  nrestart = 5,
  center = FALSE,
  scale = FALSE,
  nneg = TRUE
)
###----------------------
### d-Select number of components e.g. explaining ≥ 70%
###----------------------

# Variance explained
prop_var_addi <- peav(Y_global_scaled, nn_pca_addi$rotation, scale. = FALSE)
cum_var_addi <- cumsum(prop_var_addi)
optimal_r_nnpca_addi <- which(cum_var_addi >= 0.93)[1]
optimal_r_nnpca_addi
#optimal_r_nnpca <-optimal_r_nnpca_addi

###----------------------
### e-Extract non-negative scores and loadings
###----------------------
W_nnpca_addi <- nn_pca_addi$x[, 1:optimal_r_nnpca_addi]             # genotype × PC
H_nnpca_addi <- t(nn_pca_addi$rotation[, 1:optimal_r_nnpca_addi])   # PC × trait

#Save the contributions for later use
write.csv(H_nnpca_addi, "H_nnpca_scaled_addi.csv")
write.csv(W_nnpca_addi, "W_nnpca_scaled_addi.csv", row.names = FALSE)
###----------------------
### d-Reshape W_nmpca into y_{ilk} to get latent_array_nnpca_scaled[i, l, k] = y_{ilk}
###----------------------
#First, we extract the structure
genotypes <- unique(data_scaled$Pedigree)
env_names <- unique(data_scaled$Env)

n <- length(genotypes)
E <- length(env_names)

#Create latent array
latent_array_nnpca_scaled <- array(
  NA,
  dim = c(n, optimal_r_nnpca_addi, E),
  dimnames = list(
    genotypes,
    paste0("PC", 1:optimal_r_nnpca_addi),
    env_names
  )
)

#Fill latent_array_nnpca_scaled
#Since rows in Y_global_scaled follow the order of data_scaled, we split by environment as follows:
row_index <- 1

for (k in 1:E) {
  
  env_k <- env_names[k]
  
  rows_env <- which(data_scaled$Env == env_k)
  
  latent_array_nnpca_scaled[,,k] <- W_nnpca_addi[rows_env, ]
}

###----------------------
### e-Finally perform Finlay–Wilkinson per factor aka: MT-FW
###---------------------- 

MTFWR_results_nnpca_scaled_OLS <- list()

for (l in 1:optimal_r_nnpca_addi) {
  
  Y_factor <- latent_array_nnpca_scaled[, l, ]
  
  df_long <- data.frame(
    GEN = rep(rownames(Y_factor), times = ncol(Y_factor)),
    ENV = rep(colnames(Y_factor), each = nrow(Y_factor)),
    TRAIT = as.vector(Y_factor)
  )
  
  MTFWR_results_nnpca_scaled_OLS[[paste0("PC", l)]] <- FW::FW(
    y   = df_long$TRAIT,
    VAR = df_long$GEN,
    ENV = df_long$ENV,
    method = "OLS" # could be "Gibbs"
  )
}

#================
###----------------------
### f-Save as matrix the slope for each genotype and PCs
###---------------------- 
#Save the enrire object for later use
save(MTFWR_results_nnpca_scaled_OLS, file = "MTFWR_results_nnpca_scaled_OLS.Rdata")

pc_names <- names(MTFWR_results_nnpca_scaled_OLS)
beta_MTFWR_scaled_nnpca_OLS <- sapply(MTFWR_results_nnpca_scaled_OLS, function(res){
  res$b
})

rownames(beta_MTFWR_scaled_nnpca_OLS) <- unique(as.character(MTFWR_results_nnpca_scaled_OLS[[1]]$VAR))
write.csv(beta_MTFWR_scaled_nnpca_OLS, "beta_MTFWR_scaled_nnpca_OLS.csv")

###----------------------
### g- Comparing single tait and multi-trait
###---------------------- 

# beta_STFW_scaled_OLS: genotype × trait
# beta_MTFWR_scaled_nnpca_OLS: genotype × pc
#Pearson (linear correlation)
pearson_matrix_nnpca <- matrix(
  NA,
  nrow = ncol(beta_STFW_scaled_OLS),
  ncol = ncol(beta_MTFWR_scaled_nnpca_OLS)
)

for (t in 1:ncol(beta_STFW_scaled_OLS)) {
  for (l in 1:ncol(beta_MTFWR_scaled_nnpca_OLS)) {
    pearson_matrix_nnpca[t, l] <- cor(
      beta_STFW_scaled_OLS[, t],
      beta_MTFWR_scaled_nnpca_OLS[, l],
      method = "pearson"
    )
  }
}

rownames(pearson_matrix_nnpca) <- colnames(beta_STFW_scaled_OLS)
colnames(pearson_matrix_nnpca) <- colnames(beta_MTFWR_scaled_nnpca_OLS)

pearson_matrix_nnpca
#===============================================

#Combining all correlations such that the resulting correlation matrix contain:
#Top-left block: pc–pc plasticity correlation #Bottom-right block: Trait–trait plasticity correlation #Off-diagonal blocks: Trait–pc alignment

# Check that row order is identical
all(rownames(beta_STFW_scaled_OLS) == rownames(beta_MTFWR_scaled_nnpca_OLS))

# Combine the two slope matrices
beta_combined_scaled_nnpca_OLS <- cbind(beta_MTFWR_scaled_nnpca_OLS, beta_STFW_scaled_OLS)

#Pearson full matrix
cor_pearson_full_scaled_nnpca_OLS <- cor(beta_combined_scaled_nnpca_OLS, method = "pearson")
cor_pearson_full_scaled_nnpca_OLS
#save

write.csv(cor_pearson_full_scaled_nnpca_OLS, "cor_pearson_full_scaled_nnpca_OLS.csv")

### ---------------------------------------------------------
#Subseting trait-pcs only
#Identify pc and trait columns

pc_cols <- colnames(beta_MTFWR_scaled_nnpca_OLS)
trait_cols  <- colnames(beta_STFW_scaled_OLS)

#Extract matrices from the combined matrix

beta_MTFWR_scaled_nnpca_OLS <- beta_combined_scaled_nnpca_OLS[, pc_cols]
beta_STFW_scaled_OLS  <- beta_combined_scaled_nnpca_OLS[, trait_cols]

#Compute Pearson correlation (traits × pcs)
cor_pearson_traits_pcs_scaled_nnpca_OLS <- cor(
  beta_STFW_scaled_OLS,
  beta_MTFWR_scaled_nnpca_OLS,
  method = "pearson",
  use = "pairwise.complete.obs"
)

write.csv(cor_pearson_traits_pcs_scaled_nnpca_OLS, "cor_pearson_traits_pcs_scaled_nnpca_OLS.csv")

#==================================================================================
### ---------------------------------------------------------
### C. MTFWR Implementation on min-max scaled data with NMF
### ---------------------------------------------------------

#Using the object "data_scaled" list of matrices for all considered env
#Recalling that for each element:
#rows = genotypes
#columns = traits
#first column = Genotype ID

###----------------------
### a-Stack environments/done
###----------------------

# Identify trait columns
trait_cols <- colnames(data_scaled)[-(1:5)]

# Extract numeric matrix
Y_global_scaled <- as.matrix(data_scaled[, trait_cols])

# Ensure non-negativity (important for NMF)
if (any(Y_global_scaled < 0)) stop("Negative values detected!")

#------------------------------------------------------------------
# Choose r such that ≥70% variance explained this is the number of metatraits used to approximate the target matrix
#One could use the cophenetic correlation coefficient ( Brunet et al) as the rule of chosing the optimal r.
#REf for NMF: https://link.springer.com/article/10.1186/1471-2105-11-367
#------------------------------------------------------------------

###----------------------
### b-Select the optimal rank for factorization
###----------------------

###----------------------
###----------------------
#Uncomment if needed
#optimal_r_nmf_cophenic <- which.min(get_optimal_r_nmf$measures$rss)
#optimal_r_nmf_cophenic <- which.max(get_optimal_r_nmf$measures$cophenetic)
#optimal_r_nmf_cophenic = 3 #we got 3 for the maximum cophenic corr 
#more details at: https://cran.r-project.org/web/packages/NMF/NMF.pdf
###----------------------
###----------------------

#For consistency with nnpca, we set rank selection based on evar(explained variance) 
# optimal_r_nmf_evar <- which(get_optimal_r_nmf$measures$evar >= 0.95)[1] #depending on the data and objective, change the fraction


###----------------------
### c-Fit the final NMF
###----------------------

#optimal_r_nmf <- optimal_r_nmf_cophenic
#optimal_r_nmf <- optimal_r_nmf_evar
# nmf_final <- nmf(Y_global_scaled, rank = optimal_r_nmf, nrun = 30, seed=1234)
# 
# W_nmf_evar <- basis(nmf_final)
# H_nmf_evar <- coef(nmf_final)
# 
# #Save the contributions for later use
# write.csv(H_nmf_evar, "H_nmf_evar_scaled.csv")

###----------------------
### c-Alternative fit for NMF using the package :scNMF

#Maybe remove actual yield related trait from Y_global_scaled?
# set.seed(1234)
# nmf_final_scnmf <- scNMF::nnmf(
#   Y_global_scaled,
#   k = optimal_r_nmf,
#   rel.tol = 0.001,
#   n.threads = 0,
#   verbose = TRUE,
#   trace = 5
# )
# 
# W_nmf_scnmf <- nmf_final_scnmf$W
# H_nmf_scnmf_scaled <- nmf_final_scnmf$H
# write.csv(H_nmf_scnmf_scaled, "H_nmf_scnmf_scaled.csv")
###----------------------
##Another alternative package (What is used in the study)
library(RcppML)
#Using "mse" loss specific for Gaussian data and it is the default

nmf_final_RcppML <- RcppML::nmf(
  Y_global_scaled,
  k = 9,
  tol = 1e-05,
  maxit = 500,
  verbose = TRUE,
  L1 = .5, # L1 =c(0, 1),  # #L21 = c(0.1, 0)
  seed = 1234,
  mask_zeros = FALSE,
  diag = TRUE,
  nonneg = TRUE
)

W_nmf_RcppML <- nmf_final_RcppML$w
H_nmf_RcppML_scaled <- nmf_final_RcppML$h

colnames(H_nmf_RcppML_scaled) <- colnames(data_scaled)[-(1:5)]
rownames(H_nmf_RcppML_scaled) <- paste0("Factor", " ", 1:nrow(H_nmf_RcppML_scaled))
write.csv(H_nmf_RcppML_scaled, "H_nmf_RcppML_scaled.csv")

colnames(W_nmf_RcppML) <- paste0("Factor", " ", 1:ncol(W_nmf_RcppML))
W_nmf_RcppML <- as.data.frame(W_nmf_RcppML)
W_nmf_RcppML <- cbind(Pedigree = data_scaled$Pedigree, W_nmf_RcppML)
write.csv(W_nmf_RcppML, "W_nmf_RcppML_scaled.csv", row.names = FALSE)
###----------------------
### d-Reshape W_nmf into y_{ilk} to get latent_array_nmf_scaled[i, l, k] = y_{ilk}
###----------------------
#First, we extract the structure
optimal_r_nmf <- nrow(H_nmf_RcppML_scaled)
genotypes <- unique(data_scaled$Pedigree)
env_names <- unique(data_scaled$Env)

n <- length(genotypes)
E <- length(env_names)

#Create latent array
latent_array_nmf_scaled <- array(
  NA,
  dim = c(n, optimal_r_nmf, E),
  dimnames = list(
    genotypes,
    paste0("Factor", 1:optimal_r_nmf),
    env_names
  )
)
#Fill latent_array_nmf_scaled
#Since rows in Y_global_scaled follow the order of data_scaled, we split by environment as follows:
row_index <- 1
W_nmf_RcppML <- nmf_final_RcppML$w
colnames(W_nmf_RcppML) <- paste0("Factor", " ", 1:ncol(W_nmf_RcppML))
for (k in 1:E) {
  
  env_k <- env_names[k]
  
  rows_env <- which(data_scaled$Env == env_k)
  
  latent_array_nmf_scaled[,,k] <- W_nmf_RcppML[rows_env, ]
}

###----------------------
### e-Finally perform Finlay–Wilkinson per factor aka: MTFWR
###---------------------- 

MTFWR_results_nmf_scaled_OLS <- list()

for (l in 1:optimal_r_nmf) {
  
  Y_factor <- latent_array_nmf_scaled[, l, ]
  
  df_long <- data.frame(
    GEN = rep(rownames(Y_factor), times = ncol(Y_factor)),
    ENV = rep(colnames(Y_factor), each = nrow(Y_factor)),
    TRAIT = as.vector(Y_factor)
  )
  
  MTFWR_results_nmf_scaled_OLS[[paste0("Factor", l)]] <- FW::FW(
    y   = df_long$TRAIT,
    VAR = df_long$GEN,
    ENV = df_long$ENV,
    method = "OLS"
  )
}

###----------------------
### f-Save as matrix the slope for each genotype and factors
###---------------------- 
#Save the enrire object for later use
save(MTFWR_results_nmf_scaled_OLS, file = "MTFWR_results_nmf_scaled_OLS.Rdata")

beta_MTFWR_scaled_nmf_OLS <- sapply(MTFWR_results_nmf_scaled_OLS, function(res){
  res$b
})

rownames(beta_MTFWR_scaled_nmf_OLS) <- unique(as.character(MTFWR_results_nmf_scaled_OLS[[1]]$VAR))

#factor_names <- names(MTFWR_results_nmf_scaled_OLS)
rownames(beta_MTFWR_scaled_nmf_OLS) <- unique(as.character(MTFWR_results_nmf_scaled_OLS[[1]]$VAR))
write.csv(beta_MTFWR_scaled_nmf_OLS, "beta_MTFWR_scaled_nmf_OLS.csv")


###----------------------
### g- Comparing single tait and multi-trait
###---------------------- 

# beta_STFW_scaled_OLS: genotype × trait
# beta_MTFWR_scaled_nmf_OLS: genotype × factor
#Pearson (linear correlation)
pearson_matrix_nmf <- matrix(
  NA,
  nrow = ncol(beta_STFW_scaled_OLS),
  ncol = ncol(beta_MTFWR_scaled_nmf_OLS)
)

for (t in 1:ncol(beta_STFW_scaled_OLS)) {
  for (l in 1:ncol(beta_MTFWR_scaled_nmf_OLS)) {
    pearson_matrix_nmf[t, l] <- cor(
      beta_STFW_scaled_OLS[, t],
      beta_MTFWR_scaled_nmf_OLS[, l],
      method = "pearson"
    )
  }
}

rownames(pearson_matrix_nmf) <- colnames(beta_STFW_scaled_OLS)
colnames(pearson_matrix_nmf) <- colnames(beta_MTFWR_scaled_nmf_OLS)

pearson_matrix_nmf


#===============================================
#Combining all correlations such that the resulting correlation matrix contain:
#Top-left block: Factor–factor plasticity correlation #Bottom-right block: Trait–trait plasticity correlation #Off-diagonal blocks: Trait–factor alignment

# Check that row order is identical
all(rownames(beta_STFW_scaled_OLS) == rownames(beta_MTFWR_scaled_nmf_OLS))

# Combine the two slope matrices
beta_combined_scaled_nmf_OLS <- cbind(beta_MTFWR_scaled_nmf_OLS, beta_STFW_scaled_OLS)

#Pearson full matrix
cor_pearson_full_scaled_nmf_OLS <- cor(beta_combined_scaled_nmf_OLS, method = "pearson")
cor_pearson_full_scaled_nmf_OLS
#save

write.csv(cor_pearson_full_scaled_nmf_OLS, "cor_pearson_full_scaled_nmf_OLS.csv")

#Subseting trait-factors only
#Identify factor and trait columns
factor_cols <- colnames(beta_MTFWR_scaled_nmf_OLS)
trait_cols  <- colnames(beta_STFW_scaled_OLS)

#Extract matrices from the combined matrix

beta_MTFWR_scaled_nmf_OLS <- beta_combined_scaled_nmf_OLS[, factor_cols]
beta_STFW_scaled_OLS  <- beta_combined_scaled_nmf_OLS[, trait_cols]

#Compute Pearson correlation (traits × factors)
cor_pearson_traits_factors_scaled_nmf_OLS <- cor(
  beta_STFW_scaled_OLS,
  beta_MTFWR_scaled_nmf_OLS,
  method = "pearson",
  use = "pairwise.complete.obs"
)

write.csv(cor_pearson_traits_factors_scaled_nmf_OLS, "cor_pearson_traits_factors_scaled_nmf_OLS.csv")

#==================================================================================
#==================================================================================
#==================================================================================
#==================================================================================
#VISUALIZATION USING PEARSON AND NNPCA 
library(tidyverse)
cor_pearson_trait_pc <- cor_pearson_traits_pcs_scaled_nnpca_OLS
# transpose loadings to match trait × PC structure
loading_trait_pc <- t(H_nnpca_addi)

# ensure same ordering
loading_trait_pc <- loading_trait_pc[rownames(cor_pearson_trait_pc), colnames(cor_pearson_trait_pc)]

#Build visualization matrix
plot_df_pc_pearson <- expand.grid(
  Trait = rownames(cor_pearson_trait_pc),
  PC = colnames(cor_pearson_trait_pc)
)

plot_df_pc_pearson$correlation <- as.vector(cor_pearson_trait_pc)
plot_df_pc_pearson$loading <- as.vector(loading_trait_pc)

#Plot with ggplot
library(ggplot2)

ggplot(plot_df_pc_pearson, aes(x = PC, y = Trait)) +
  
  geom_tile(aes(fill = correlation), color = "white") +
  
  geom_point(aes(size = abs(loading)), color = "black", alpha = 0.7) +
  
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = bquote(rho(beta[ST], beta[PC]))
  ) +
  
  scale_size_continuous(
    name = "Trait loading"
  ) +
  
  # theme_minimal(base_size = 14) +
  
  labs(
    #title = "Relationship between STFW-OLS plasticity and MTFW PC plasticity",
    #x = " ",
    x = "Metatraits (MTFW)",
    y = "Traits (STFW)"
  ) +
  
  theme(plot.title.position = 'plot', 
        text = element_text(size = 11),
        plot.title  = element_text(hjust = 0.5),
        axis.text.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"))

#================================================
#================================================
#PLot loading in the x axis and cor betaST, betaMT in y axis
#direct diagnostic plot showing, for each Trait–PC pair as per Zo's suggestion
#x=loadingtrait,PC
#y=cor(βtrait,βPC)
#This could served as sanity checks for multi-trait Finlay–Wilkinson models because it tests:
#if or not: β_PC≈∑_j H_PC,j βj
#The question is: are traits with high loading also showing strong correlation with the PC plasticity?
#We will include statistical test for pair-wise correlation

### ---------------------------------------------------------
### 1. Prepare the data, compute correlations and p-values.
### Because each correlation uses the same genotypes, we aim to compute Pearson/Kendall correlation
### then, tests for each Trait–PC pair and correct for multiple testing
### We have: 
### beta_STFW_scaled_OLS: genotype × traits
### beta_MTFWR_scaled_nnpca_OLS: genotype × PCs
### ---------------------------------------------------------
traits <- colnames(beta_STFW_scaled_OLS)
pcs    <- colnames(beta_MTFWR_scaled_nnpca_OLS)

cor_mat_pearson_nnpca  <- matrix(NA, nrow = length(traits), ncol = length(pcs))
pval_mat_pearson_nnpca <- matrix(NA, nrow = length(traits), ncol = length(pcs))

rownames(cor_mat_pearson_nnpca)  <- traits
colnames(cor_mat_pearson_nnpca)  <- pcs
rownames(pval_mat_pearson_nnpca) <- traits
colnames(pval_mat_pearson_nnpca) <- pcs

for (t in traits) {
  for (p in pcs) {
    
    test <- cor.test(
      beta_STFW_scaled_OLS[,t],
      beta_MTFWR_scaled_nnpca_OLS[,p],
      method = "pearson"
    )
    
    cor_mat_pearson_nnpca[t,p]  <- test$estimate
    pval_mat_pearson_nnpca[t,p] <- test$p.value
  }
}

### ---------------------------------------------------------
### 2. Correct for multiple testing
### ---------------------------------------------------------
pval_adj_pearson_nnpca <- matrix(
  p.adjust(pval_mat_pearson_nnpca, method = "BH"),
  nrow = nrow(pval_mat_pearson_nnpca)
)

rownames(pval_adj_pearson_nnpca) <- rownames(pval_mat_pearson_nnpca)
colnames(pval_adj_pearson_nnpca) <- colnames(pval_mat_pearson_nnpca)

### ---------------------------------------------------------
### 3. Build dataframe for plotting
### ---------------------------------------------------------
library(tidyverse)

loading_trait_pc <- t(H_nnpca_addi)

loading_trait_pc <- loading_trait_pc[rownames(cor_mat_pearson_nnpca),colnames(cor_mat_pearson_nnpca)]

plot_df_pc_pearson <- expand.grid(
  Trait = rownames(cor_mat_pearson_nnpca),
  PC = colnames(cor_mat_pearson_nnpca)
)

plot_df_pc_pearson$correlation <- as.vector(cor_mat_pearson_nnpca)
plot_df_pc_pearson$loading     <- as.vector(loading_trait_pc)
plot_df_pc_pearson$pvalue      <- as.vector(pval_adj_pearson_nnpca)
### ---------------------------------------------------------
### 4. Select and highlight traits with significant correlations and loadings
### ---------------------------------------------------------
#This is for significant correlation only
# label_df <- plot_df_pc_pearson %>%
#   filter(pvalue < 0.05)
# label_df <- plot_df_pc_pearson %>%
#   filter(pvalue < 0.05 & abs(loading) > 0.2 & abs(correlation) >= 0.5)

label_df <- plot_df_pc_pearson %>%
  filter(pvalue < 0.05 & abs(correlation) >= 0.5)
### ---------------------------------------------------------
### 5. Visualization with significance highlighted
### ---------------------------------------------------------
library(ggplot2)
library(ggrepel)

### ---------------------------------------------------------
### Using the "plot_df_pc_pearson" (Trait–PC pairs) object
### ---------------------------------------------------------
### ### a. compute alignment score
### ---------------------------------------------------------
library(dplyr)
library(ggplot2)

### ---------------------------------------------------------
### Compute alignment score and p-value
### ---------------------------------------------------------

alignment_scores_pearson_nnpca <- plot_df_pc_pearson %>%
  group_by(PC) %>%
  summarise(
    alignment = cor(loading, correlation, use = "complete.obs"),
    pvalue = cor.test(loading, correlation)$p.value,
    .groups = "drop"
  )

### ---------------------------------------------------------
### FDR correction
### ---------------------------------------------------------

alignment_scores_pearson_nnpca$FDR <- p.adjust(
  alignment_scores_pearson_nnpca$pvalue,
  method = "BH"
)

### ---------------------------------------------------------
### Define significance
### require BOTH statistical significance AND effect size
### ---------------------------------------------------------

alignment_scores_pearson_nnpca <- alignment_scores_pearson_nnpca %>%
  mutate(
    signif = case_when(
      FDR <= 0.001 & abs(alignment) >= 0.5 ~ "***",
      FDR <= 0.01  & abs(alignment) >= 0.5 ~ "**",
      FDR <= 0.05  & abs(alignment) >= 0.5 ~ "*",
      TRUE ~ ""
    )
  )

### ---------------------------------------------------------
### Label position
### ---------------------------------------------------------

alignment_scores_pearson_nnpca <- alignment_scores_pearson_nnpca %>%
  mutate(label_pos = alignment + 0.05)

### ---------------------------------------------------------
### 5. Save results
### ---------------------------------------------------------

write.csv(
  alignment_scores_pearson_nnpca,
  "alignment_scores_nnpca_pearson_OLS.csv",
  row.names = FALSE
)

### ---------------------------------------------------------
### Plot
### ---------------------------------------------------------

ggplot(alignment_scores_pearson_nnpca,
       aes(x = PC, y = alignment)) +
  
  geom_col(fill = "steelblue") +
  
  geom_text(
    aes(label = round(alignment, 2)),
    vjust = -0.3,
    size = 4
  ) +
  
  geom_text(
    aes(y = label_pos, label = signif),
    size = 6,
    color = "red"
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  
  expand_limits(
    y = max(alignment_scores_pearson_nnpca$label_pos) + 0.05
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  
  labs(
    #y = "Plasticity alignment score",
    y = bquote("Plasticity alignment score" ~ (rho)), 
    x = " "
    #x = "Principal components"
  )

#===============================================================================================
##Alignment between trait importance and plasticity 

#Ranking still based on loading magnitude
plot_df_pc_pearson <- plot_df_pc_pearson %>%
  group_by(PC) %>%
  arrange(desc(abs(loading)), .by_group = TRUE) %>%
  mutate(
    Rank = row_number(),
    Loading_abs = abs(loading)
  ) %>%
  ungroup()
### ---------------------------------------------------------
#Keep signed correlations
threshold <- 0.5

plot_df_pc_pearson <- plot_df_pc_pearson %>%
  mutate(
    strong = abs(correlation) >= threshold
  )
### ---------------------------------------------------------
#plot (preserve sign clearly)
ggplot(plot_df_pc_pearson, aes(x = Rank, y = correlation)) +
  
  # baseline
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  # all points
  geom_point(color = "grey70", alpha = 0.5) +
  
  # highlight strong correlations
  geom_point(
    data = subset(plot_df_pc_pearson, strong),
    aes(color = correlation),
    size = 2
  ) +
  
  # smooth trend
  geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
  
  facet_wrap(~PC, scales = "free_x") +
  
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  
  theme_minimal(base_size = 14) +
  
  labs(
    x = "Trait ranking by loading",
    y =  bquote(rho(beta[ST], beta[PC])),
    color = "Correlation",
    #title = "Alignment between trait importance and plasticity (signed)",
    #subtitle = "Red = positive alignment, Blue = opposite response"
  )
### ---------------------------------------------------------
#Remove the grid
ggplot(plot_df_pc_pearson, aes(x = Rank, y = correlation)) +
  
  # baseline
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  # all points
  geom_point(color = "grey70", alpha = 0.5) +
  
  # highlight strong correlations
  geom_point(
    data = subset(plot_df_pc_pearson, strong),
    aes(color = correlation),
    size = 2
  ) +
  
  # smooth trend
  geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
  
  facet_wrap(~PC, scales = "free_x") +
  
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    panel.grid = element_blank(),        # remove all grid
    panel.border = element_rect(         # optional: add clean border
      color = "black", fill = NA, linewidth = 0.5
    ),
    axis.line = element_line(color = "black")
  ) +
  
  labs(
    x = "Trait ranking by loading",
    y =  bquote(rho(beta[ST], beta[PC])),
    color = "Correlation",
    #title = "Alignment between trait importance and plasticity (signed)",
    #subtitle = "Red = positive alignment, Blue = opposite response"
  )

#==================================================================================
#==================================================================================
#==================================================================================
#==================================================================================
#VISUALIZATION USING PEARSON AND NMF 
library(tidyverse)
cor_pearson_trait_factor <- cor_pearson_traits_factors_scaled_nmf_OLS
# transpose loadings to match trait × Factor structure
loading_trait_factor <- t(H_nmf_RcppML_scaled)

# ensure same ordering
#loading_trait_factor <- loading_trait_factor[rownames(cor_pearson_trait_factor), colnames(cor_pearson_trait_factor)]

#Build visualization matrix
plot_df_factor_pearson <- expand.grid(
  Trait = rownames(cor_pearson_trait_factor),
  Factor = colnames(cor_pearson_trait_factor)
)

plot_df_factor_pearson$correlation <- as.vector(cor_pearson_trait_factor)
plot_df_factor_pearson$loading <- as.vector(loading_trait_factor)

#Plot with ggplot
library(ggplot2)

ggplot(plot_df_factor_pearson, aes(x = Factor, y = Trait)) +
  
  geom_tile(aes(fill = correlation), color = "white") +
  
  geom_point(aes(size = abs(loading)), color = "black", alpha = 0.7) +
  
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = bquote(rho(beta[ST], beta[Factor]))
  ) +
  
  scale_size_continuous(
    name = "Trait loading"
  ) +
  
  # theme_minimal(base_size = 14) +
  
  labs(
    #title = "Relationship between STFW-OLS plasticity and MTFW Factor plasticity",
    x = "Metatraits (MTFW)",
    y = "Traits (STFW)"
  ) +
  
  theme(plot.title.position = 'plot', 
        text = element_text(size = 11),
        plot.title  = element_text(hjust = 0.5),
        axis.text.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"))

### ---------------------------------------------------------
### 1. Prepare the data, compute correlations and p-values.
### Because each correlation uses the same genotypes, we aim to compute Pearson/Kendall correlation
### then, tests for each Trait–Factor pair and correct for multiple testing
### We have: 
### beta_STFW_scaled_OLS: genotype × traits
### beta_MTFWR_scaled_nmf_OLS: genotype × Factors
### ---------------------------------------------------------
traits <- colnames(beta_STFW_scaled_OLS)
factors    <- colnames(beta_MTFWR_scaled_nmf_OLS)

cor_mat_pearson_nmf  <- matrix(NA, nrow = length(traits), ncol = length(factors))
pval_mat_pearson_nmf <- matrix(NA, nrow = length(traits), ncol = length(factors))

rownames(cor_mat_pearson_nmf)  <- traits
colnames(cor_mat_pearson_nmf)  <- factors
rownames(pval_mat_pearson_nmf) <- traits
colnames(pval_mat_pearson_nmf) <- factors

for (t in traits) {
  for (p in factors) {
    
    test <- cor.test(
      beta_STFW_scaled_OLS[,t],
      beta_MTFWR_scaled_nmf_OLS[,p],
      method = "pearson"
    )
    
    cor_mat_pearson_nmf[t,p]  <- test$estimate
    pval_mat_pearson_nmf[t,p] <- test$p.value
  }
}

### ---------------------------------------------------------
### 2. Correct for multiple testing
### ---------------------------------------------------------
pval_adj_pearson_nmf <- matrix(
  p.adjust(pval_mat_pearson_nmf, method = "BH"),
  nrow = nrow(pval_mat_pearson_nmf)
)

rownames(pval_adj_pearson_nmf) <- rownames(pval_mat_pearson_nmf)
colnames(pval_adj_pearson_nmf) <- colnames(pval_mat_pearson_nmf)

### ---------------------------------------------------------
### 3. Build dataframe for plotting
### ---------------------------------------------------------
library(tidyverse)

loading_trait_factor <- t(H_nmf_RcppML_scaled)

#loading_trait_factor <- loading_trait_factor[rownames(cor_mat_pearson_nmf),colnames(cor_mat_pearson_nmf)]

plot_df_factor_pearson <- expand.grid(
  Trait = rownames(cor_mat_pearson_nmf),
  Factor = colnames(cor_mat_pearson_nmf)
)

plot_df_factor_pearson$correlation <- as.vector(cor_mat_pearson_nmf)
plot_df_factor_pearson$loading     <- as.vector(loading_trait_factor)
plot_df_factor_pearson$pvalue      <- as.vector(pval_adj_pearson_nmf)
### ---------------------------------------------------------
### 4. Select and highlight traits with significant correlations and loadings
### ---------------------------------------------------------
#This is for significant correlation only
# label_df <- plot_df_factor_pearson %>%
#   filter(pvalue < 0.05)
# label_df <- plot_df_factor_pearson %>%
#   filter(pvalue < 0.05 & abs(loading) > 0.2 & abs(correlation) >= 0.5)

label_df <- plot_df_factor_pearson %>%
  filter(pvalue < 0.05 & abs(correlation) >= 0.5)
### ---------------------------------------------------------
### 5. Visualization with significance highlighted
### ---------------------------------------------------------
library(ggplot2)
library(ggrepel)
### ---------------------------------------------------------
### 6. Show significance with symbols
### ---------------------------------------------------------
plot_df_factor_pearson$signif <- cut(
  plot_df_factor_pearson$pvalue,
  breaks = c(-Inf,0.001,0.01,0.05,Inf),
  labels = c("***","**","*","")
)

#Label them
aes(label = paste0(Trait, " ", signif))
### ---------------------------------------------------------
### ---------------------------------------------------------
### Using the "plot_df_factor_pearson" (Trait–Factor pairs) object
### ---------------------------------------------------------
### ### a. compute alignment score
### ---------------------------------------------------------
library(dplyr)
library(ggplot2)

### ---------------------------------------------------------
### Compute alignment score and p-value
### ---------------------------------------------------------

alignment_scores_pearson_nmf <- plot_df_factor_pearson %>%
  group_by(Factor) %>%
  summarise(
    alignment = cor(loading, correlation, use = "complete.obs"),
    pvalue = cor.test(loading, correlation)$p.value,
    .groups = "drop"
  )

### ---------------------------------------------------------
### FDR correction
### ---------------------------------------------------------

alignment_scores_pearson_nmf$FDR <- p.adjust(
  alignment_scores_pearson_nmf$pvalue,
  method = "BH"
)

### ---------------------------------------------------------
### Define significance
### require BOTH statistical significance AND effect size
### ---------------------------------------------------------

alignment_scores_pearson_nmf <- alignment_scores_pearson_nmf %>%
  mutate(
    signif = case_when(
      FDR <= 0.001 & abs(alignment) >= 0.5 ~ "***",
      FDR <= 0.01  & abs(alignment) >= 0.5 ~ "**",
      FDR <= 0.05  & abs(alignment) >= 0.5 ~ "*",
      TRUE ~ ""
    )
  )

### ---------------------------------------------------------
### Label position
### ---------------------------------------------------------

alignment_scores_pearson_nmf <- alignment_scores_pearson_nmf %>%
  mutate(label_pos = alignment + 0.05)

### ---------------------------------------------------------
### 5. Save results
### ---------------------------------------------------------

write.csv(
  alignment_scores_pearson_nmf,
  "alignment_scores_nmf_pearson_OLS.csv",
  row.names = FALSE
)

### ---------------------------------------------------------
### Plot
### ---------------------------------------------------------

ggplot(alignment_scores_pearson_nmf,
       aes(x = Factor, y = alignment)) +
  
  geom_col(fill = "steelblue") +
  
  geom_text(
    aes(label = round(alignment, 2)),
    vjust = -0.3,
    size = 4
  ) +
  
  geom_text(
    aes(y = label_pos, label = signif),
    size = 6,
    color = "red"
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  
  expand_limits(
    y = max(alignment_scores_pearson_nmf$label_pos) + 0.05
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  
  labs(
    y = bquote("Plasticity alignment score" ~ (rho)),
    x = " "
  )

#===============================================================================================
#===============================================================================================

##Alignment between trait importance and plasticity 

#Ranking still based on loading magnitude
plot_df_factor_pearson <- plot_df_factor_pearson %>%
  group_by(Factor) %>%
  arrange(desc(abs(loading)), .by_group = TRUE) %>%
  mutate(
    Rank = row_number(),
    Loading_abs = abs(loading)
  ) %>%
  ungroup()
### ---------------------------------------------------------
#Keep signed correlations
threshold <- 0.5

plot_df_factor_pearson <- plot_df_factor_pearson %>%
  mutate(
    strong = abs(correlation) >= threshold
  )
### ---------------------------------------------------------
#plot (preserve sign clearly)
ggplot(plot_df_factor_pearson, aes(x = Rank, y = correlation)) +
  
  # baseline
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  # all points
  geom_point(color = "grey70", alpha = 0.5) +
  
  # highlight strong correlations
  geom_point(
    data = subset(plot_df_factor_pearson, strong),
    aes(color = correlation),
    size = 2
  ) +
  
  # smooth trend
  geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
  
  facet_wrap(~Factor, scales = "free_x") +
  
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  
  theme_minimal(base_size = 14) +
  
  labs(
    x = "Trait ranking by loading",
    y =  bquote(rho(beta[ST], beta[Factor])),
    color = "Correlation"
  )
### ---------------------------------------------------------

# ===============================================
# Trait–Factor Concordance Analysis 
# ===============================================

library(dplyr)
library(readr)
library(Hmisc)

# --------------------------
# LOAD DATA
# --------------------------

# Trait × Factor correlation matrix
cor_trait_fac <- read.csv(
  "cor_pearson_traits_factors_scaled_nmf_OLS.csv",
  row.names = 1,
  check.names = FALSE
)

# Phenotype data
pheno <- read.csv(
  "Maize_Pheno_Data_4Env.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Remove metadata columns
trait_cols <- colnames(pheno)[-(1:5)]
pheno <- pheno[, trait_cols]

# --------------------------
# PARAMETERS
# --------------------------
threshold_tf <- 0.5   # trait-factor
threshold_tt <- 0.5   # trait-trait
alpha <- 0.05         # significance

# --------------------------
# STEP 1: Select traits per factor
# --------------------------
traits_per_factor <- lapply(colnames(cor_trait_fac), function(fac) {
  traits <- rownames(cor_trait_fac)[abs(cor_trait_fac[, fac]) > threshold_tf]
  return(traits)
})

names(traits_per_factor) <- colnames(cor_trait_fac)

# --------------------------
# STEP 2: Trait–Trait correlation analysis
# --------------------------
analyze_trait_correlations <- function(traits, pheno_data) {
  
  traits <- intersect(traits, colnames(pheno_data))
  
  if(length(traits) < 2) return(NULL)
  
  df <- pheno_data[, traits, drop = FALSE]
  df <- df[complete.cases(df), ]
  
  if(nrow(df) < 5) return(NULL)
  
  cor_res <- Hmisc::rcorr(as.matrix(df), type = "pearson")
  
  cor_mat <- cor_res$r
  p_mat <- cor_res$P
  
  results <- data.frame()
  
  for(i in 1:(ncol(df)-1)){
    for(j in (i+1):ncol(df)){
      
      results <- rbind(results, data.frame(
        Trait1 = colnames(df)[i],
        Trait2 = colnames(df)[j],
        Correlation = cor_mat[i, j],
        P_value = p_mat[i, j]
      ))
    }
  }
  
  # FDR correction
  results$FDR <- p.adjust(results$P_value, method = "fdr")
  
  # --------------------------
  # STRICT FILTER (FIX)
  # --------------------------
  results <- results %>%
    filter(
      abs(Correlation) > threshold_tt,
      FDR < alpha
    )
  
  return(results)
}

# --------------------------
# STEP 3: Loop over factors
# --------------------------
results_list_fac <- list()

for(fac in names(traits_per_factor)){
  
  traits <- traits_per_factor[[fac]]
  
  if(length(traits) < 2) next
  
  res <- analyze_trait_correlations(traits, pheno)
  
  if(is.null(res) || nrow(res) == 0) next
  
  res$Factor <- fac
  
  results_list_fac[[fac]] <- res
}

# --------------------------
# STEP 4: Combine results
# --------------------------
final_results_fac <- bind_rows(results_list_fac)

# --------------------------
# SAVE OUTPUT
# --------------------------
write.csv(
  final_results_fac,
  "Trait_Factor_Concordance_Analysis_Maize146.csv",
  row.names = FALSE
)

# ===============================================
# Trait–PC Concordance Analysis 
# ===============================================
# --------------------------
# LOAD DATA
# --------------------------

# Trait × PC correlation matrix
cor_trait_pc <- read.csv(
  "cor_pearson_traits_pcs_scaled_nnpca_OLS.csv",
  row.names = 1,
  check.names = FALSE
)

# Phenotype data
pheno <- read.csv(
  "Maize_Pheno_Data_4Env.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Remove metadata columns
trait_cols <- colnames(pheno)[-(1:5)]
pheno <- pheno[, trait_cols]

# --------------------------
# PARAMETERS
# --------------------------
threshold_tp <- 0.5   # trait-pc
threshold_tt <- 0.5   # trait-trait
alpha <- 0.05         # significance

# --------------------------
# STEP 1: Select traits per pc
# --------------------------
traits_per_pc <- lapply(colnames(cor_trait_pc), function(pc) {
  traits <- rownames(cor_trait_pc)[abs(cor_trait_pc[, pc]) > threshold_tp]
  return(traits)
})

names(traits_per_pc) <- colnames(cor_trait_pc)

# --------------------------
# STEP 2: Trait–Trait correlation analysis
# --------------------------
analyze_trait_correlations <- function(traits, pheno_data) {
  
  traits <- intersect(traits, colnames(pheno_data))
  
  if(length(traits) < 2) return(NULL)
  
  df <- pheno_data[, traits, drop = FALSE]
  df <- df[complete.cases(df), ]
  
  if(nrow(df) < 5) return(NULL)
  
  cor_res <- Hmisc::rcorr(as.matrix(df), type = "pearson")
  
  cor_mat <- cor_res$r
  p_mat <- cor_res$P
  
  results <- data.frame()
  
  for(i in 1:(ncol(df)-1)){
    for(j in (i+1):ncol(df)){
      
      results <- rbind(results, data.frame(
        Trait1 = colnames(df)[i],
        Trait2 = colnames(df)[j],
        Correlation = cor_mat[i, j],
        P_value = p_mat[i, j]
      ))
    }
  }
  
  # FDR correction
  results$FDR <- p.adjust(results$P_value, method = "fdr")
  
  # --------------------------
  # STRICT FILTER 
  # --------------------------
  results <- results %>%
    filter(
      abs(Correlation) > threshold_tt,
      FDR < alpha
    )
  
  return(results)
}

# --------------------------
# STEP 3: Loop over pcs
# --------------------------
results_list_pc <- list()

for(pc in names(traits_per_pc)){
  
  traits <- traits_per_pc[[pc]]
  
  if(length(traits) < 2) next
  
  res <- analyze_trait_correlations(traits, pheno)
  
  if(is.null(res) || nrow(res) == 0) next
  
  res$PC <- pc
  
  results_list_pc[[pc]] <- res
}

# --------------------------
# STEP 4: Combine results
# --------------------------
final_results_pc <- bind_rows(results_list_pc)

# --------------------------
# SAVE OUTPUT
# --------------------------
write.csv(
  final_results_pc,
  "Trait_PC_Concordance_Analysis_Maize146.csv",
  row.names = FALSE
)
### ---------------------------------------------------------
### ---------------------------------------------------------

