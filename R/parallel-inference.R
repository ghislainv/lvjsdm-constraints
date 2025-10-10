require(jSDM)
require(doParallel)
require(foreach)

parallel_inference <- function(Y, X, nchains, burnin, mcmc, thin,
                               starting_values, n_latent, seed) {
  ## Make a cluster for parallel MCMCs
  ncores <- nchains ## One core for each MCMC chains
  clust <- makeCluster(ncores)
  registerDoParallel(clust)

  ## Parallel inference
  models <- foreach(i=1:nchains) %dopar% {
    mod <- jSDM::jSDM_binomial_probit(
      # Iteration
      burnin=burnin,
      mcmc=mcmc,
      thin=thin,
      # Response variable
      presence_data=Y,
      # Explanatory variables
      site_formula=~x1+x2,
      site_data=X,
      n_latent=n_latent,
      site_effect="random",
      # Starting values
      alpha_start=starting_values$alpha[i],
      beta_start=starting_values$beta[i],
      lambda_start=getElement(starting_values$lambda, i),
      W_start=starting_values$W[i],
      V_alpha=starting_values$V_alpha[i],
      # Priors
      shape_Valpha=0.5,
      rate_Valpha=0.0005,
      mu_beta=0, V_beta=1,
      mu_lambda=0, V_lambda=1,
      seed=seed, verbose=1)
    return(mod)
  }

  # Stop the cluster
  stopCluster(clust)

  # Return the model
  return(models)
}
