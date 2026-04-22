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
library(DHARMa)
library(glue)

# Source files
source("R/sim-data.R")
source("R/parallel-inference.R")
source("R/convergence.R")

# Model names
model_names <- c("nosort", "sort", "Z", "dharma", "PIT", "qres")
n_models <- 6

# =======================================
# Simulate data
# =======================================

out_dir <- here("outputs")
seed <- 1234
data <- sim_data_LVJSDM(
  n_species=30,
  n_sites=100,
  n_p=3,
  n_q=3,
  V_alpha_target=0.1,
  neg_val_diag=-0.1,
  seed=seed,
  out_dir=out_dir)
X <- data$X
Y <- data$Y
W <- data$W
beta_target <- data$beta_target
lambda_target <- data$lambda_target
n_sites <- nrow(X)
n_q <- ncol(X)
n_species <- ncol(Y)
ofile <- file.path(out_dir, "lambda_target.rda")
save(lambda_target, file=ofile)

# =======================================
# MCMC parameters
# =======================================

burnin <- 10000
mcmc <- 10000
thin <- 10
ngibbs <- burnin + mcmc
nchains <- 8

# =======================================
# Priors
# =======================================

# Prior variance for lambdas
V_lambda <- 1

# =======================================
# Starting values
# =======================================

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

# ===========================================
# Model 1 no sorting --> convergence problems
# ===========================================

# mod_nosort
ofile <- file.path(out_dir, "mod_nosort.rda")
if (!file.exists(ofile)) {
  mod_nosort <- parallel_inference(
    Y, X,
    nchains=nchains,
    burnin=burnin, mcmc=mcmc, thin=thin,
    n_latent=n_q, V_lambda=V_lambda,
    starting_values=starting_values, seed=seed)
  save(mod_nosort, file=ofile)
} else {
  load(ofile)
}

# Rhat
Rhat_nosort <- compute_rhat(mod_nosort) |>
  rename(maxRhat_nosort=maxRhat) |>
  write_csv(file.path(out_dir, "rhat_nosort.csv"))

# Sum of median VCV and Cor
vcv_cor_median_nosort <- compute_vcv_cor_median(
  mod_nosort, id_target_species=c(1, 2, 3)) |>
  rename(sum_nosort=sum) |>
  write_csv(file.path(out_dir, "vcv_cor_med_nosort.csv"))

# Plot traces
mcmc_list <- get_mcmc_list_lambdas(
  mod_nosort,
  re="(sp_1\\.lambda_1|sp_4\\.lambda_1)"
)
png(here(out_dir, "mcmc_nosort.png"))
plot(mcmc_list)
dev.off()

# =======================================
# Model 2 sorting species manually
# =======================================

# Sorting species
sp_max <- c(4, 5, 6)
sp_max_sort <- sp_max
Y <- read.csv(
  file=file.path(out_dir, "Y.csv"),
  header=TRUE, row.names=1)
Y_sort <- Y
Y_sort[, c(1, 2, 3)] <- Y[, sp_max]
Y_sort[, sp_max] <- Y[, c(1, 2, 3)]

# mod_sort
ofile <- file.path(out_dir, "mod_sort.rda")
if (!file.exists(ofile)) {
  mod_sort <- parallel_inference(
    Y_sort, X,
    nchains=nchains,
    burnin=burnin, mcmc=mcmc, thin=thin,
    n_latent=n_q, V_lambda=V_lambda,
    starting_values=starting_values, seed=seed)
  save(mod_sort, file=ofile)
} else {
  load(ofile)
}

# Rhat
Rhat_sort <- compute_rhat(mod_sort) |>
  rename(maxRhat_sort=maxRhat) |>
  write_csv(file.path(out_dir, "rhat_sort.csv"))

# Sum of median VCV and Cor
id_tsp <- sp_max_sort
vcv_cor_median_sort <- compute_vcv_cor_median(
  mod_sort, id_target_species=id_tsp) |>
  rename(sum_sort=sum) |>
  write_csv(file.path(out_dir, "vcv_cor_med_sort.csv"))

# Plot traces
mcmc_list <- get_mcmc_list_lambdas(
  mod_sort,
  re="(sp_1\\.lambda_1|sp_4\\.lambda_1)"
)
png(file.path(out_dir, "mcmc_sort.png"))
plot(mcmc_list)
dev.off()

## # Identify incorrect chains when starting lambdas != 0
## unlist(lapply(mcmc_list, function(x) {mean(x[,1])}))
## ## --> chains 5 and 8
## # Initial values
## W_start
## lambda_start[c(4, 8)]

# =======================================
# Model with no latent variables
# =======================================

# Starting values
starting_values_nolv <- starting_values[c("alpha", "V_alpha", "beta")]

# mod_nolv (no latent variables)
mod_nolv <- parallel_inference(
  Y, X,
  nchains=1,
  burnin=burnin, mcmc=mcmc, thin=thin,
  n_latent=0,
  starting_values=starting_values_nolv,
  seed=seed)
mod_nolv <- mod_nolv[[1]]

# =======================================
# Model 3 sorting species automatically
# using PCA on residuals (feasable with jSDM, not HSMC)
# =======================================

# Residuals from latent Z
Z_latent <- mod_nolv$Z_latent  # ! Predictive posterior means
aXbeta <- matrix(0, nrow=n_sites, ncol=n_species)
for (i in 1:n_sites) {
  X_i <- X[i, ]
  alpha_i <- mod_nolv$mcmc.alpha[, i]
  for (j in 1:n_species) {
    beta_j <- mod_nolv$mcmc.sp[[j]]
    aXbeta[i, j] <- mean(alpha_i + X_i %*% t(beta_j))
  }
}
residuals <- Z_latent - aXbeta

# ACP on residuals
pca <- prcomp(residuals,  scale=FALSE)
sp_max_PC1 <- which.max(abs(pca$rotation[, "PC1"]))
sp_max_PC2 <- which.max(abs(pca$rotation[, "PC2"]))
sp_max_PC3 <- which.max(abs(pca$rotation[, "PC3"]))
sp_max <- c(sp_max_PC1, sp_max_PC2, sp_max_PC3)
sp_max_Z <- sp_max

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
ggsave(file.path(out_dir, "loadings_pca_lv_Z.png"), p)
corr <- loadings |>
  group_by(axis) |>
  summarize(corr=round(cor(pca, lv), 2)) |>
  ungroup() |>
  write_csv(file.path(out_dir, "loadings_corr_Z.csv"))

# Sorting species
Y <- read.csv(
  file=file.path(out_dir, "Y.csv"),
  header=TRUE, row.names=1)
Y_sort_pca <- Y
Y_sort_pca[, c(1, 2, 3)] <- Y[, sp_max]
Y_sort_pca[, sp_max] <- Y[, c(1, 2, 3)]

# mod_Z
ofile <- file.path(out_dir, "mod_Z.rda")
if (!file.exists(ofile)) {
  mod_Z <- parallel_inference(
    Y_sort_pca, X,
    nchains=nchains,
    burnin=burnin, mcmc=mcmc, thin=thin,
    n_latent=n_q, V_lambda=V_lambda,
    starting_values=starting_values, seed=seed)
  save(mod_Z, file=ofile)
} else {
  load(ofile)
}

# Rhat
Rhat_Z <- compute_rhat(mod_Z) |>
  rename(maxRhat_Z=maxRhat) |>
  write_csv(file.path(out_dir, "rhat_Z.csv"))

# Sum of median VCV and Cor
id_tsp <- sp_max_Z
vcv_cor_median_Z <- compute_vcv_cor_median(
  mod_Z, id_target_species=id_tsp) |>
  rename(sum_Z=sum) |>
  write_csv(file.path(out_dir, "vcv_cor_med_Z.csv"))

# Plot traces
mcmc_list <- get_mcmc_list_lambdas(
  mod_Z,
  re="(sp_1\\.lambda_1|sp_4\\.lambda_1)"
)
png(file.path(out_dir, "mcmc_Z.png"))
plot(mcmc_list)
dev.off()

# ===========================================
# Model 4 sorting species automatically
# using PCA on quantile residuals from DHARMa
# ===========================================

# Simulate data for each observation
# from model with no latent variables.
n_sim <- mcmc / thin
sim_data <- matrix(NA, nrow=n_sites * n_species, ncol=n_sim)
theta_sim <- array(NA, dim=c(n_sites, n_species, n_sim))
# Prepare beta_sim
beta_sim <- array(NA, dim=c(n_q, n_species, n_sim))
for (j in 1:n_species) {
  # Transpose mcmc to get dimensions c(n_q, n_sim)
  beta_sim[, j, ] <- t(mod_nolv$mcmc.sp[[j]])
}
# Set seed and simulate observations
set.seed(seed)
for (s in 1:n_sim) {
  alpha_sim <- mod_nolv$mcmc.alpha[s, ]
  e_sim <- matrix(rnorm(n_sites * n_species, 0, 1), n_sites, n_species)
  Z_sim <- alpha_sim + X %*% beta_sim[, , s] + e_sim
  # pnorm is inverse-probit
  theta_sim[, , s] <- matrix(pnorm(c(Z_sim)), nrow=n_sites)
  # Presence-absence matrix Y
  Y_sim <- matrix(0, n_sites, n_species)
  Y_sim[Z_sim > 0] <- 1
  sim_data[, s] <- c(Y_sim)
}
## # Check that theta_median ~= sim_mean
## theta_median <- apply(theta_sim, c(1, 2), median)
## sim_mean <- apply(sim_data, 1, mean)
## plot(c(theta_median), sim_mean)
## # Check residuals
## dharma <- createDHARMa(
##   simulatedResponse=sim_data,
##   observedResponse=c(as.matrix(Y)),
##   fittedPredictedResponse=c(theta_median),
##   integerResponse=TRUE, seed=seed)
## plot(dharma)

# Get residuals using package DHARMa
qres <- DHARMa::getQuantile(
  simulations=sim_data,
  observed=c(as.matrix(Y)),
  integerResponse=TRUE,
  method="traditional")

# Scaled residuals as matrix nsites * nspecies
qresiduals <- matrix(qres, nrow=n_sites)

# ACP on residuals
pca <- prcomp(qresiduals,  scale=FALSE)
sp_max_PC1 <- which.max(abs(pca$rotation[, "PC1"]))
sp_max_PC2 <- which.max(abs(pca$rotation[, "PC2"]))
sp_max_PC3 <- which.max(abs(pca$rotation[, "PC3"]))
sp_max <- c(sp_max_PC1, sp_max_PC2, sp_max_PC3)
sp_max_dharma <- sp_max

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
ggsave(file.path(out_dir, "loadings_pca_lv_dharma.png"), p)
corr <- loadings |>
  group_by(axis) |>
  summarize(corr=round(cor(pca, lv), 2)) |>
  ungroup() |>
  write_csv(file.path(out_dir, "loadings_corr_dharma.csv"))

# Sorting species
Y <- read.csv(
  file=file.path(out_dir, "Y.csv"),
  header=TRUE, row.names=1)
Y_sort_pca <- Y
Y_sort_pca[, c(1, 2, 3)] <- Y[, sp_max]
Y_sort_pca[, sp_max] <- Y[, c(1, 2, 3)]

# mod_dharma
ofile <- file.path(out_dir, "mod_dharma.rda")
if (!file.exists(ofile)) {
  mod_dharma <- parallel_inference(
    Y_sort_pca, X,
    nchains=nchains,
    burnin=burnin, mcmc=mcmc, thin=thin,
    n_latent=n_q, V_lambda=V_lambda,
    starting_values=starting_values, seed=seed)
  save(mod_dharma, file=ofile)
} else {
  load(ofile)
}
    
# Rhat
Rhat_dharma <- compute_rhat(mod_dharma) |>
  rename(maxRhat_dharma=maxRhat) |>
  write_csv(file.path(out_dir, "rhat_dharma.csv"))

# Sum of median VCV and Cor
id_tsp <- sp_max_dharma
vcv_cor_median_dharma <- compute_vcv_cor_median(
  mod_dharma, id_target_species=id_tsp) |>
  rename(sum_dharma=sum) |>
  write_csv(file.path(out_dir, "vcv_cor_med_dharma.csv"))

# Plot traces
mcmc_list <- get_mcmc_list_lambdas(
  mod_dharma,
  re="(sp_1\\.lambda_1|sp_4\\.lambda_1)"
)
png(file.path(out_dir, "mcmc_dharma.png"))
plot(mcmc_list)
dev.off()

# ===========================================
# Model 5 sorting species automatically
# using PCA on PIT residuals
# ===========================================

## PIT: probability integral transform.
## For a reference on PIT residuals, see Warton et al. 2017:
## https://doi.org/10.1371/journal.pone.0181790

## For implementing PIT algorithm, see:
## https://github.com/florianhartig/DHARMa/issues/168#issuecomment-623155563

# Simulations and observations
simulations <- sim_data
observed <- c(as.matrix(Y))

# Compute PIT residuals
min <- apply(simulations < observed, 1, sum) / n_sim
max <- apply(simulations <= observed, 1, sum) / n_sim
PIT_residuals <- runif(n_sites * n_species, min, max)

# Residuals as matrix nsites * nspecies
qresiduals <- matrix(PIT_residuals, nrow=n_sites)

# ACP on residuals
pca <- prcomp(qresiduals,  scale=FALSE)
sp_max_PC1 <- which.max(abs(pca$rotation[, "PC1"]))
sp_max_PC2 <- which.max(abs(pca$rotation[, "PC2"]))
sp_max_PC3 <- which.max(abs(pca$rotation[, "PC3"]))
sp_max <- c(sp_max_PC1, sp_max_PC2, sp_max_PC3)
sp_max_PIT <- sp_max

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
ggsave(file.path(out_dir, "loadings_pca_lv_PIT.png"), p)
corr <- loadings |>
  group_by(axis) |>
  summarize(corr=round(cor(pca, lv), 2)) |>
  ungroup() |>
  write_csv(file.path(out_dir, "loadings_corr_PIT.csv"))

# Sorting species
Y <- read.csv(
  file=file.path(out_dir, "Y.csv"),
  header=TRUE, row.names=1)
Y_sort_pca <- Y
Y_sort_pca[, c(1, 2, 3)] <- Y[, sp_max]
Y_sort_pca[, sp_max] <- Y[, c(1, 2, 3)]

# mod_PIT
ofile <- file.path(out_dir, "mod_PIT.rda")
if (!file.exists(ofile)) {
  mod_PIT <- parallel_inference(
    Y_sort_pca, X,
    nchains=nchains,
    burnin=burnin, mcmc=mcmc, thin=thin,
    n_latent=n_q, V_lambda=V_lambda,
    starting_values=starting_values, seed=seed)
  save(mod_PIT, file=ofile)
} else {
  load(ofile)
}
  
# Rhat
Rhat_PIT <- compute_rhat(mod_PIT) |>
  rename(maxRhat_PIT=maxRhat) |>
  write_csv(file.path(out_dir, "rhat_PIT.csv"))

# Sum of median VCV and Cor
id_tsp <- sp_max_PIT
vcv_cor_median_PIT <- compute_vcv_cor_median(
  mod_PIT, id_target_species=id_tsp) |>
  rename(sum_PIT=sum) |>
  write_csv(file.path(out_dir, "vcv_cor_med_PIT.csv"))

# Plot traces
mcmc_list <- get_mcmc_list_lambdas(
  mod_PIT,
  re="(sp_1\\.lambda_1|sp_4\\.lambda_1)"
)
png(file.path(out_dir, "mcmc_PIT.png"))
plot(mcmc_list)
dev.off()

# ===========================================
# Model 6 sorting species automatically
# using n_sim PCA on quantile residuals computed manually
# ===========================================

# Observations
y <- c(as.matrix(Y))

# Predicted probabilities
p_pred <- matrix(c(theta_sim), ncol=n_sim)

# 1. Calculate the bounds of the interval [a, b]
# pbinom(k, size, prob) is the CDF of the binomial distribution
a <- pbinom(y - 1, size = 1, prob = p_pred)
b <- pbinom(y, size = 1, prob = p_pred)

# 2. Random sampling (Jittering)
# Draw a uniform value u within the interval [a, b]
# This handles the discreteness of the binomial distribution
u <- matrix(runif(length(y) * n_sim, min = a, max = b), ncol=n_sim)

# 3. Transformation via the normal quantile function (Probit)
# qnorm transforms probabilities into normal distribution quantiles
qres <- qnorm(u)

# 4. PCA on residuals and species with highest coordinates
sp_max_mat <- matrix(NA, nrow=3, ncol=n_sim)
for (i in 1:n_sim) {
  # Residuals as matrix nsites * nspecies
  qresiduals <- matrix(qres[, i], ncol=n_species)
  # ACP on residuals
  pca <- prcomp(qresiduals,  scale=FALSE)
  sp_max_PC1 <- which.max(abs(pca$rotation[, "PC1"]))
  sp_max_PC2 <- which.max(abs(pca$rotation[, "PC2"]))
  sp_max_PC3 <- which.max(abs(pca$rotation[, "PC3"]))
  sp_max_mat[, i] <- c(sp_max_PC1, sp_max_PC2, sp_max_PC3)
}
# Function to get the mode of a vector 
get_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
# Get species with highest coordinates on average
sp_max <- apply(sp_max_mat, 1, get_mode)
sp_max_qres <- sp_max

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
ggsave(file.path(out_dir, "loadings_pca_lv_qres.png"), p)
corr <- loadings |>
  group_by(axis) |>
  summarize(corr=round(cor(pca, lv), 2)) |>
  ungroup() |>
  write_csv(file.path(out_dir, "loadings_corr_qres.csv"))

# Sorting species
Y <- read.csv(
  file=file.path(out_dir, "Y.csv"),
  header=TRUE, row.names=1)
Y_sort_pca <- Y
Y_sort_pca[, c(1, 2, 3)] <- Y[, sp_max]
Y_sort_pca[, sp_max] <- Y[, c(1, 2, 3)]

# mod_qres
ofile <- file.path(out_dir, "mod_qres.rda")
if (!file.exists(ofile)) {
  mod_qres <- parallel_inference(
    Y_sort_pca, X,
    nchains=nchains,
    burnin=burnin, mcmc=mcmc, thin=thin,
    n_latent=n_q, V_lambda=V_lambda,
    starting_values=starting_values, seed=seed)
  save(mod_qres, file=ofile)
} else {
  load(ofile)
}
 
# Rhat
Rhat_qres <- compute_rhat(mod_qres) |>
  rename(maxRhat_qres=maxRhat) |>
  write_csv(file.path(out_dir, "rhat_qres.csv"))

# Sum of median VCV and Cor
id_tsp <- sp_max_qres
vcv_cor_median_qres <- compute_vcv_cor_median(
  mod_qres, id_target_species=id_tsp) |>
  rename(sum_qres=sum) |>
  write_csv(file.path(out_dir, "vcv_cor_med_qres.csv"))

# Plot traces
mcmc_list <- get_mcmc_list_lambdas(
  mod_qres,
  re="(sp_1\\.lambda_1|sp_4\\.lambda_1)"
)
png(file.path(out_dir, "mcmc_qres.png"))
plot(mcmc_list)
dev.off()

# =======================================
# Approach comparison
# =======================================

# Rhat
# ---

# Load Rhat tables if necessary
for (i in 1:n_models) {
  mod <- model_names[i]
  if (!glue("Rhat_{mod}") %in% ls()) {
    ifile <- file.path(out_dir, glue("rhat_{mod}.csv"))
    identity_transformer(glue("Rhat_{mod} <- read_csv(ifile)"))
  }
}

Rhat_df <- Rhat_nosort |>
  left_join(Rhat_sort, by="Variable") |>
  left_join(Rhat_Z, by="Variable") |>
  left_join(Rhat_dharma, by="Variable") |>
  left_join(Rhat_PIT, by="Variable") |>
  left_join(Rhat_qres, by="Variable") |>
  write_csv(file.path(out_dir, "Rhat_model_comparison.csv"))

# VCV, Cor
# --------

# Load vcv/cor tables if necessary
for (i in 1:n_models) {
  mod <- model_names[i]
  if (!glue("vcv_cor_median_{mod}") %in% ls()) {
    ifile <- file.path(out_dir, glue("vcv_cor_med_{mod}.csv"))
    identity_transformer(glue("vcv_cor_median_{mod} <- read_csv(ifile)"))
  }
}

# From lambda_target
vcv_cor_target <- compute_vcov_cor(t(lambda_target))
vcv_target <- vcv_cor_target$vcv
cor_target <- vcv_cor_target$cor
id_tsp <- c(1, 2, 3)
sum_var_target <- sum(diag(vcv_target[id_tsp, id_tsp]))
sum_cov_target <- sum(vcv_target[id_tsp, ]) - sum_var_target
sum_cor_target <- sum(cor_target[id_tsp, ])
var_names <- c("V_tsp", "absCov_tsp", "absCorr_tsp")
sum_values <- round(c(sum_var_target, sum_cov_target,
                      sum_cor_target), 2)
vcv_cor_target <- data.frame(Variable=var_names, sum_target=sum_values)

vcv_cor_median_df <- vcv_cor_median_nosort |>
  left_join(vcv_cor_median_sort, by="Variable") |>
  left_join(vcv_cor_median_Z, by="Variable") |>
  left_join(vcv_cor_median_dharma, by="Variable") |>
  left_join(vcv_cor_median_PIT, by="Variable") |>
  left_join(vcv_cor_median_qres, by="Variable") |>
  left_join(vcv_cor_target, by="Variable") |>
  write_csv(file.path(out_dir, "vcv_cor_median_model_comparison.csv"))

# =======================================
# Plot correlations estimates vs. target
# =======================================

# ------------------------
# For lambdas
# ------------------------

load(file=file.path(out_dir, "lambda_target.rda"))
lambda_123 <- lambda_target[, c(1, 2, 3)] 

# 1. mod_nosort
p_nosort <- plot_corr_comp(mod_nosort, lambda_target,
                           id_species=c(1, 2, 3), "mean") +
  ggtitle("Model nosort")

# 2-6 models
for (i in 2:n_models) {
  mod <- model_names[i]
  sp_max <- identity_transformer(glue("sp_max_{mod}"))
  lambda_target_sort <- lambda_target
  lambda_target_sort[, c(1, 2, 3)] <- lambda_target[, sp_max]
  lambda_target_sort[, sp_max] <- lambda_123
  model <- identity_transformer(glue("mod_{mod}"))
  p <- plot_corr_comp(model, lambda_target_sort,
                      id_species=c(1, 2, 3), "mean") +
    ggtitle(glue("Model {mod}"))
  identity_transformer(glue("p_{mod} <- p"))
}

library(patchwork)
p_all <- p_nosort + p_sort + p_Z + p_dharma + p_PIT + p_qres + plot_layout(ncol=3)
p_all

# End
