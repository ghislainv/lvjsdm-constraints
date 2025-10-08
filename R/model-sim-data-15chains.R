## ==============================================================================
## author          :Ghislain Vieilledent
## email           :ghislain.vieilledent@cirad.fr
## web             :https://ecology.ghislainv.fr
## license         :GPLv3
## ==============================================================================

## Libraries
library(jSDM)
library(doParallel)
library(foreach)
library(dplyr)
library(here)

## Source files
source("R/sim-data.R")
source("R/parallel-inference.R")
source("R/convergence.R")

## =======================================
## Simulate data
## =======================================

out_dir <- here("outputs")
data <- sim_data_LVJSDM(n_species=30,
                        n_sites=100,
                        n_p=3,  # Number of explanatory variables (including intercept)
                        n_q=3,  # Number of latent variables
                        V_alpha_target=0.1,  # Variance of random site effects
                        neg_val_diag=-0.1, # Small negative values on the diagonal
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

## =======================================
## Starting values
## =======================================

nchains <- 15
seed <- 1234
set.seed(seed)
alpha_start <- rnorm(nchains, 0, 1)
beta_start <- rnorm(nchains, 0, 1)
lambda_start <- rnorm(nchains, 0, 1)
W_start <- rnorm(nchains, 0, 1)
V_alpha_start <- runif(nchains, 1, 2)
starting_values <- list(
  alpha=alpha_start, beta=beta_start, lambda=lambda_start,
  W=W_start, V_alpha=V_alpha_start)

## =======================================
## MCMC parameters
## =======================================

burnin <- 5000
mcmc <- 5000
thin <- 5
ngibbs <- burnin + mcmc
nchains <- 5

## =======================================
## Model 1
## =======================================

## mod_1
mod_1 <- parallel_inference(
  Y, X,
  nchains=nchains,
  burnin=burnin, mcmc=mcmc, thin=thin,
  starting_values=starting_values, seed=seed)

## Rhat
Rhat_1 <- compute_rhat(mod_1)

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
Y <- Y_sort

## mod_2
mod_2 <- parallel_inference(
  Y, X,
  nchains=nchains,
  burnin=burnin, mcmc=mcmc, thin=thin,
  starting_values=starting_values, seed=seed)

## Rhat
Rhat_2 <- compute_rhat(mod_2)

# =======================================
# Model 3 sorting species automatically
# =======================================



# =======================================
# Approach comparison
# =======================================

Rhat_df <- Rhat_1 |>
  left_join(Rhat_2, by="Variable") |>
  rename(maxRhat_no=maxRhat.x, maxRhat_to=maxRhat.y)

## End
