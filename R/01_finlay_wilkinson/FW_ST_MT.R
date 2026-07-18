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
#------------------------------------------------------------------

###----------------------
### b-Select the optimal rank for factorization
###----------------------

###----------------------
###----------------------
#Uncomment if needed
#optimal_r_nmf_cophenic <- which.min(get_optimal_r_nmf$measures$rss)
#optimal_r_nmf_cophenic <- which.max(get_optimal_r_nmf$measures$cophenetic)
#more details at: https://cran.r-project.org/web/packages/NMF/NMF.pdf
###----------------------
###----------------------

#For consistency with nnpca, we set rank selection based on evar(explained variance) 
# optimal_r_nmf_evar <- which(get_optimal_r_nmf$measures$evar >= 0.95)[1] #depending on the data and objective, change the fraction


###----------------------
### c-Fit the final NMF
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
