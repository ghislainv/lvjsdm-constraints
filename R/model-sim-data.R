# ===============================================
# author          :Ghislain Vieilledent
# email           :ghislain.vieilledent@cirad.fr
# web             :https://ecology.ghislainv.fr
# license         :GPLv3
# ===============================================

# Libraries
library(jSDM)
library(doParallel)
library(foreach)
library(dplyr)
library(here)
library(ggplot2)
library(readr)

# Source files
source("R/sim-data.R")
source("R/parallel-inference.R")
source("R/convergence.R")

# =======================================
# Simulate data
# =======================================

out_dir <- here("outputs")
data <- sim_data_LVJSDM(
  n_species=30,
  n_sites=100,
  n_p=3,
  n_q=3,
  V_alpha_target=0.1,
  neg_val_diag=-0.1,
  seed=1234,
  out_dir=out_dir)
X <- data$X
Y <- data$Y
W <- data$W
beta_target <- data$beta_target
lambda_target <- data$lambda_target
n_sites <- nrow(X)
n_q <- ncol(X)
n_species <- ncol(Y)

# =======================================
# MCMC parameters
# =======================================

burnin <- 10000
mcmc <- 10000
thin <- 10
ngibbs <- burnin + mcmc
nchains <- 8

# =======================================
# Starting values
# =======================================

seed <- 1234
set.seed(seed)
alpha_start <- rnorm(nchains, 0, 1)
beta_start <- rnorm(nchains, 0, 1)
## # lambda_start as list of matrices
## lambda_start <- list()
## lambda_start_mat <- matrix(0, n_q, n_species)
## axis_imp <- seq(n_q, 1, -1)  ## Decreasing importance of the latent axis
## for (i in 1:nchains) {
##   mat0 <- matrix(runif(n_species * n_q, -0.2, 0.2), nrow=n_q, ncol=n_species)
##   mat <- mat0 * axis_imp
##   lambda_start_mat[upper.tri(mat, diag=TRUE)] <- mat[upper.tri(mat, diag=TRUE)]
##   diag(lambda_start_mat) <- axis_imp * runif(1, 0.5, 1.5)
##   lambda_start[[i]] <- lambda_start_mat
## }
lambda_start <- rep(0, nchains) # Set to zero for convergence
W_start <- rnorm(nchains, 0, 1)
V_alpha_start <- runif(nchains, 0, 1)
starting_values <- list(
  alpha=alpha_start, beta=beta_start, lambda=lambda_start,
  W=W_start, V_alpha=V_alpha_start)

# =======================================
# Model 1
# =======================================

# mod_1
mod_1 <- parallel_inference(
  Y, X,
  nchains=nchains,
  burnin=burnin, mcmc=mcmc, thin=thin,
  n_latent=n_q,
  starting_values=starting_values, seed=seed)

# Rhat
Rhat_1 <- compute_rhat(mod_1)

# Plot traces
mcmc_list <- mcmc_lambdas(
  mod_1,
  re="(sp_1\\.lambda_1|sp_4\\.lambda_1)"
)
png(here(out_dir, "mcmc_nosort.png"))
plot(mcmc_list)
dev.off()

# =======================================
# Model 2 sorting species manually
# =======================================

# Sorting species
Y <- read.csv(file=file.path(out_dir, "Y.csv"), header=TRUE, row.names=1)
Y_sort <- Y
Y_sort[, 1] <- Y[, 4]
Y_sort[, 2] <- Y[, 5]
Y_sort[, 3] <- Y[, 6]
Y_sort[, 4] <- Y[, 1]
Y_sort[, 5] <- Y[, 2]
Y_sort[, 6] <- Y[, 3]

# mod_2
mod_2 <- parallel_inference(
  Y_sort, X,
  nchains=nchains,
  burnin=burnin, mcmc=mcmc, thin=thin,
  n_latent=n_q,
  starting_values=starting_values, seed=seed)

# Rhat
Rhat_2 <- compute_rhat(mod_2)

# Plot traces
mcmc_list <- mcmc_lambdas(
  mod_2,
  re="(sp_1\\.lambda_1|sp_4\\.lambda_1)"
)
png(here(out_dir, "mcmc_sort.png"))
plot(mcmc_list)
dev.off()

## # Identify incorrect chains
## unlist(lapply(mcmc_list, function(x) {mean(x[,1])}))
## ## --> chains 5 and 8
## # Initial values
## W_start
## lambda_start[c(4, 8)]

# =======================================
# Model 3 sorting species automatically
# using PCA on residuals
# =======================================

# Starting values
starting_values_3 <- starting_values[c("alpha", "V_alpha", "beta")]

# mod_3
mod_3 <- parallel_inference(
  Y, X,
  nchains=1,
  burnin=burnin, mcmc=mcmc, thin=thin,
  n_latent=0,
  starting_values=starting_values_3, seed=seed)
mod_3 <- mod_3[[1]]

# Residuals
Z_latent <- mod_3$Z_latent  # ! Predictive posterior means
aXbeta <- matrix(0, nrow=n_sites, ncol=n_species)
for (i in 1:n_sites) {
  X_i <- X[i, ]
  alpha_i <- mod_3$mcmc.alpha[, i]
  for (j in 1:n_species) {
    beta_j <- mod_3$mcmc.sp[[j]]
    aXbeta[i, j] <- mean(alpha_i + X_i %*% t(beta_j))
  }
}
residuals <- Z_latent - aXbeta

# ACP
pca <- prcomp(residuals,  scale=FALSE)
sp_max_PC1 <- which.max(abs(pca$rotation[, "PC1"]))
sp_max_PC2 <- which.max(abs(pca$rotation[, "PC2"]))
sp_max_PC3 <- which.max(abs(pca$rotation[, "PC3"]))

# Correlation between PCA and LV loadings
loadings <- data.frame(
  pca=c(pca$rotation[, c("PC1", "PC2", "PC3")]),
  lv=c(t(lambda_target)),
  axis=rep(c("Axis 1", "Axis 2", "Axis 3"), each=n_species))
p <- loadings |>
  ggplot(aes(lv, pca, col=axis)) +
  geom_point() +
  facet_grid(rows=vars(axis)) +
  xlab("lambda targets") +
  ylab("PCA loadings") +
  theme_bw() + guides(colour="none")
ggsave(here(out_dir, "loadings_pca_lv.png"), p)
corr <- loadings |>
  group_by(axis) |>
  summarize(corr=round(cor(pca, lv), 2)) |>
  ungroup() |>
  write_csv(here(out_dir, "loadings_corr.csv"))

# Sorting species
Y <- read.csv(
  file=file.path(out_dir, "Y.csv"),
  header=TRUE, row.names=1)
Y_sort_pca <- Y
Y_sort_pca[, 1] <- Y[, sp_max_PC1]
Y_sort_pca[, 2] <- Y[, sp_max_PC2]
Y_sort_pca[, 3] <- Y[, sp_max_PC3]
Y_sort_pca[, sp_max_PC1] <- Y[, 1]
Y_sort_pca[, sp_max_PC2] <- Y[, 2]
Y_sort_pca[, sp_max_PC3] <- Y[, 3]

# mod_3
mod_3 <- parallel_inference(
  Y_sort_pca, X,
  nchains=nchains,
  burnin=burnin, mcmc=mcmc, thin=thin,
  n_latent=n_q,
  starting_values=starting_values, seed=seed)

# Rhat
Rhat_3 <- compute_rhat(mod_3)

# Plot traces
mcmc_list <- mcmc_lambdas(
  mod_2,
  re="(sp_1\\.lambda_1|sp_4\\.lambda_1)"
)
png(here(out_dir, "mcmc_sort_pca.png"))
plot(mcmc_list)
dev.off()

# =======================================
# Model 4 sorting species automatically
# using PCA on quantile residuals
# =======================================

## To be done

# =======================================
# Approach comparison
# =======================================

Rhat_df <- Rhat_1 |>
  left_join(Rhat_2, by="Variable") |>
  rename(maxRhat_no=maxRhat.x, maxRhat_to=maxRhat.y)

# End
