require(coda)
require(ggplot2)
require(dplyr)

##' @title Transform an array into a MCMC
##' @param x An array object
##' @return An mcmc object
##' @author Ghislain Vieilledent
arr2mcmc <- function(x) {
  return(mcmc(as.data.frame(x),
              start=burnin+1 , end=ngibbs, thin=thin))
}

##' @title Get MCMC for species covariances and correlations
##' @param model jSDM model output
##' @return List of MCMCs for for species covariances and
##'   correlations. List are separated between constrained and
##'   unconstrained species.
##' @author Ghislain Vieilledent
get_vcv_cor_mcmc <- function(mod) {
  # Variables
  burnin <- mod$model_spec$burnin
  mcmc <- mod$model_spec$mcmc
  thin <- mod$model_spec$thin
  ngibbs <- burnin + mcmc
  n_species <- length(mod$mcmc.sp)
  n_samp <- mcmc / thin
  n_interactions <- (n_species^2 - n_species) / 2 + n_species
  n_q <- mod$model_spec$n_latent
  n_inter_const_sp <- (n_species * n_q) - ((n_q * n_q - n_q) / 2)
  par_names <- colnames(mod$mcmc.sp[[1]])
  # VCV and Cor
  vcv_mat <- matrix(NA, n_samp, n_interactions)
  cor_mat <- matrix(NA, n_samp, n_interactions)
  for (i in 1:n_samp) {
    lv_coefs <- t(sapply(mod$mcmc.sp, "[", i, grep("lambda", par_names)))
    vcv <- lv_coefs %*% t(lv_coefs)
    vcv_mat[i, ] <- vcv[lower.tri(vcv, diag=TRUE)]
    cor <- cov2cor(vcv)
    cor_mat[i, ] <- cor[lower.tri(cor, diag=TRUE)]
  }
  # Matrices
  vcv_mat_const_sp <- vcv_mat[, 1:n_inter_const_sp]
  vcv_mat_unconst_sp <- vcv_mat[, (n_inter_const_sp + 1):n_interactions]
  cor_mat_const_sp <- cor_mat[, 1:n_inter_const_sp]
  cor_mat_unconst_sp <- cor_mat[, (n_inter_const_sp + 1):n_interactions]
  # Data-frames
  vcv_df_const_sp <- as.data.frame(vcv_mat_const_sp)
  vcv_df_unconst_sp <- as.data.frame(vcv_mat_unconst_sp)
  cor_df_const_sp <- as.data.frame(cor_mat_const_sp)
  cor_df_unconst_sp <- as.data.frame(cor_mat_unconst_sp)  
  # MCMC
  vcv_mcmc_const_sp <- mcmc(vcv_df_const_sp, burnin+1, ngibbs, thin)
  vcv_mcmc_unconst_sp <- mcmc(vcv_df_unconst_sp, burnin+1, ngibbs, thin)
  cor_mcmc_const_sp <- mcmc(cor_df_const_sp, burnin+1, ngibbs, thin)
  cor_mcmc_unconst_sp <- mcmc(cor_df_unconst_sp, burnin+1, ngibbs, thin)
  # Result
  return(list(vcv_mcmc_const_sp=vcv_mcmc_const_sp,
              vcv_mcmc_unconst_sp=vcv_mcmc_unconst_sp,
              cor_mcmc_const_sp=cor_mcmc_const_sp,
              cor_mcmc_unconst_sp=cor_mcmc_unconst_sp))
}

##' @title Compute maximum Rhat
##' @param model_output jSDM model output with parallel chains
##' @return A data-frame with maximum Rhat value per parameter type
##' @author Ghislain Vieilledent
compute_rhat <- function(model_output) {

  ## Parameters for coda object
  burnin <- model_output[[1]]$model_spec$burnin
  ngibbs <- burnin + model_output[[1]]$model_spec$mcmc
  thin <-  model_output[[1]]$model_spec$thin

  # MCMC lists
  mcmc_list_alpha <- mcmc.list(lapply(lapply(model_output,"[[","mcmc.alpha"), arr2mcmc))
  mcmc_list_V_alpha <- mcmc.list(lapply(lapply(model_output,"[[","mcmc.V_alpha"), arr2mcmc))
  mcmc_list_sp <- mcmc.list(lapply(lapply(model_output,"[[","mcmc.sp"), arr2mcmc))
  mcmc_list_lv <- mcmc.list(lapply(lapply(model_output,"[[","mcmc.latent"), arr2mcmc))
  mcmc_list_lambda <- mcmc.list(
    lapply(mcmc_list_sp[, grep("lambda", colnames(mcmc_list_sp[[1]]), value=TRUE)], arr2mcmc))
  mcmc_list_deviance <- mcmc.list(lapply(lapply(model_output,"[[","mcmc.Deviance"), arr2mcmc))

  # Get MCMC for species covariances and correlations
  # separating constrained and unconstrained species
  vcv_corr <- lapply(model_output, get_vcv_cor_mcmc)
  vcv_mcmc_const_sp <- lapply(vcv_corr, "[[", 1)
  vcv_mcmc_unconst_sp <- lapply(vcv_corr, "[[", 2)
  cor_mcmc_const_sp <- lapply(vcv_corr, "[[", 3)
  cor_mcmc_unconst_sp <- lapply(vcv_corr, "[[", 4)
  
  # psrf gelman indice 
  psrf_alpha <- max(gelman.diag(
    mcmc_list_alpha,
    multivariate=FALSE)$psrf[,2])
  psrf_V_alpha <- gelman.diag(
    mcmc_list_V_alpha)$psrf[,2]
  psrf_beta <- max(gelman.diag(
    mcmc_list_sp[, grep("beta", colnames(mcmc_list_sp[[1]]))],
    multivariate=FALSE)$psrf[,2])
  psrf_lambda <- max(gelman.diag(
    mcmc_list_lambda, multivariate=FALSE)$psrf[,2], na.rm=TRUE)
  psrf_lv <- max(gelman.diag(
    mcmc_list_lv, multivariate=FALSE)$psrf[,2 ])
  psrf_vcv_const_sp <- max(gelman.diag(
    vcv_mcmc_const_sp, multivariate=FALSE)$psrf[,2], na.rm=TRUE)
  psrf_vcv_unconst_sp <- max(gelman.diag(
    vcv_mcmc_unconst_sp, multivariate=FALSE)$psrf[,2], na.rm=TRUE)
  psrf_cor_const_sp <- max(gelman.diag(
    cor_mcmc_const_sp, multivariate=FALSE)$psrf[,2], na.rm=TRUE)
  psrf_cor_unconst_sp <- max(gelman.diag(
    cor_mcmc_unconst_sp, multivariate=FALSE)$psrf[,2], na.rm=TRUE)

  # Combining max Rhat
  max_rhat_others <- max(psrf_alpha, psrf_V_alpha, psrf_beta)
  max_rhat_vcv_corr_csp <- max(psrf_vcv_const_sp, psrf_cor_const_sp)
  max_rhat_vcv_corr_ucsp <- max(psrf_vcv_unconst_sp, psrf_cor_unconst_sp)
  max_rhat <- round(c(max_rhat_others, psrf_lambda, psrf_lv,
                max_rhat_vcv_corr_csp, max_rhat_vcv_corr_ucsp), 2)
  var_names <- c("others", "loadings", "lv", "vcv_cor_csp", "vcv_cor_ucsp")
  Rhat <- data.frame(Variable=var_names, maxRhat=max_rhat)

  # Return result
  return(Rhat)
}

##' @title Compute vcv and corr medians for constrained species 
##' @param model_output jSDM model output with parallel chains
##' @param id_target_species Vector of id for target species. Must be c(1, 2, 3) when unsorted.
##' @return A data-frame with the vcv and corr medians for constrained species 
##' @author Ghislain Vieilledent
compute_vcv_cor_median <- function(model_output, id_target_species) {

  ## Parameters for coda object
  burnin <- model_output[[1]]$model_spec$burnin
  ngibbs <- burnin + model_output[[1]]$model_spec$mcmc
  thin <-  model_output[[1]]$model_spec$thin

  # MCMC lists
  mcmc_list_sp <- mcmc.list(lapply(lapply(model_output,"[[","mcmc.sp"), arr2mcmc))
  mcmc_list_lambda <- mcmc.list(
    lapply(mcmc_list_sp[, grep("lambda", colnames(mcmc_list_sp[[1]]), value=TRUE)], arr2mcmc))

  # Get MCMC for species covariances and correlations
  # separating constrained and unconstrained species
  vcv_corr <- lapply(model_output, get_vcv_cor_mcmc)
  vcv_mcmc_const_sp <- lapply(vcv_corr, "[[", 1)
  vcv_mcmc_unconst_sp <- lapply(vcv_corr, "[[", 2)
  cor_mcmc_const_sp <- lapply(vcv_corr, "[[", 3)
  cor_mcmc_unconst_sp <- lapply(vcv_corr, "[[", 4)

  # Combine chains
  mcmc_vcv_csp <- dplyr::bind_rows(lapply(vcv_mcmc_const_sp, data.frame))
  mcmc_vcv_ucsp <- dplyr::bind_rows(lapply(vcv_mcmc_unconst_sp, data.frame))
  mcmc_cor_csp <- dplyr::bind_rows(lapply(cor_mcmc_const_sp, data.frame))
  mcmc_cor_ucsp <- dplyr::bind_rows(lapply(cor_mcmc_unconst_sp, data.frame))

  # Compute median of absolute values
  median_vcv_csp <- apply(apply(mcmc_vcv_csp, 2, abs), 2, median)
  median_vcv_ucsp <- apply(apply(mcmc_vcv_ucsp, 2, abs), 2, median)
  median_cor_csp <- apply(apply(mcmc_cor_csp, 2, abs), 2, median)
  median_cor_ucsp <- apply(apply(mcmc_cor_ucsp, 2, abs), 2, median)

  # Lower triangular matrices
  n_species <- length(model_output[[1]]$mcmc.sp)
  ltri_vcv <- ltri_cor <- matrix(NA, n_species, n_species)
  ltri_vcv[lower.tri(ltri_vcv, diag=TRUE)] <- c(median_vcv_csp, median_vcv_ucsp)
  ltri_cor[lower.tri(ltri_cor, diag=TRUE)] <- c(median_cor_csp, median_cor_ucsp)

  # Full matrices
  full_vcv <- ltri_vcv
  full_vcv[upper.tri(full_vcv, diag=TRUE)] <- t(ltri_vcv)[upper.tri(full_vcv, diag=TRUE)]
  full_cor <- ltri_cor
  full_cor[upper.tri(full_cor, diag=TRUE)] <- t(ltri_cor)[upper.tri(full_cor, diag=TRUE)]

  # Sum of variances, abs(covariances), and abs(correlations) for target species
  sum_var_target <- sum(diag(full_vcv[id_target_species, id_target_species]))
  sum_cov_target <- sum(full_vcv[id_target_species, ]) - sum_var_target
  sum_cor_target <- sum(full_cor[id_target_species, ])

  # Make data frame
  var_names <- c("V_tsp", "absCov_tsp", "absCorr_tsp")
  sum_values <- round(c(sum_var_target, sum_cov_target, sum_cor_target), 2)
  vcv_cor_median <- data.frame(Variable=var_names, sum=sum_values)

  # Return result
  return(vcv_cor_median)
}

##' @title Plot Rhat
##' @param Rhat dataframe of Rhat values
##' @return ggplot object
##' @author Ghislain Vieilledent
plot_rhat <- function(Rhat) {
  p <- ggplot(Rhat, aes(x=Variable, y=Rhat)) + 
    ggtitle("Averages of Rhat obtained for each type of parameter") +
    theme(plot.title = element_text(hjust=0.5, size=13)) +
    geom_bar(fill="skyblue", stat = "identity") +
    geom_text(aes(label=round(Rhat, 3)), vjust=0, hjust=-0.1) +
    geom_hline(yintercept=1, color='red') +
    ylim(0, max(Rhat$Rhat) + 0.2) +
    coord_flip()
  return(p)
}

##' @title Return an mcmc list for a selection of lambdas
##' @param model_output List of model output obtained from parallel
##'   computing
##' @param re Regular expression to select lambdas. Default to
##'   "lambda" for all lambdas. Could also be
##'   "(sp_1\\.lambda_1|sp_4\\.lambda_1)" for example.
##' @return
##' @author Ghislain Vieilledent
get_mcmc_list_lambdas <- function(model_output, re="lambda") {
  # Parameters for coda object
  burnin <- model_output[[1]]$model_spec$burnin
  ngibbs <- burnin + model_output[[1]]$model_spec$mcmc
  thin <-  model_output[[1]]$model_spec$thin
  # MCMC lists
  mcmc_list_sp <- mcmc.list(
    lapply(lapply(model_output, "[[", "mcmc.sp"), arr2mcmc)
  )
  wcols <- grep(re, colnames(mcmc_list_sp[[1]]), value=TRUE)
  mcmc_list_selection <- mcmc_list_sp[, wcols]
  return(mcmc_list_selection)
}

##' @title Compute VCV and correlations matrices from lambdas
##' @param lambdas Species loadings. Matrix of dimensions n_species * n_latent.
##' @return VCV and correlations matrix
##' @author Ghislain Vieilledent
compute_vcov_cor <- function(lambdas) {
  vcv <- lambdas %*% t(lambdas)
  cor <- cov2cor(vcv)
  return(list(vcv=vcv, cor=cor))
}

##' @title Tranform an mcmc list into a dataframe
##' @param mcmc_list
##' @return dataframe
##' @author Ghislain Vieilledent
mcmc_list2df <- function(mcmc_list) {
  list_df <- lapply(mcmc_list, as.data.frame)
  df <- dplyr::bind_rows(list_df)
  return(df)
}

##' @title Comparing estimated and target species correlations
##' @param model_output Model output with several chains
##' @param lambda_target Lambda target values. Matrix with dimensions
##'   n_q * n_species
##' @param id_species Id of constrained species in the unsorted model
##' @param estimate "mean" to compute mean or something else for median
##' @return ggplot with correlation value
##' @author Ghislain Vieilledent
plot_corr_comp <- function(model_output, lambda_target,
                           id_species, estimate="mean") {
  mcmc_list_lambdas <- get_mcmc_list_lambdas(model_output)
  df_lambdas <- mcmc_list2df(mcmc_list_lambdas)
  if (estimate == "mean") {
    lambda_est <- apply(df_lambdas, 2, mean)
  } else {
    lambda_est <- apply(df_lambdas, 2, median)
  }
  lambda_t <- c(lambda_target)
  lambda_names <- names(lambda_est)
  cor_value <- round(cor(lambda_t, lambda_est), 2)
  slope <- ifelse(cor_value >= 0, 1, -1)
  id_sp <- id_species
  lambda_par <- vector()
  for (i in id_sp) {
    par <- paste0("sp_", i, paste0(".lambda_", 1:3))
    lambda_par <- c(lambda_par, par)
  }
  csp <- lambda_names %in% lambda_par
  df_plot_lambda <- data.frame(
    lambda_est,
    lambda_t,
    lambda_names,
    csp,
    row.names=NULL)
  # Plot
  xrng <- range(lambda_t)
  yrng <- range(lambda_est)
  xymax <- ceiling(max(abs(xrng), abs(yrng)))
  p <- tibble(df_plot_lambda) |>
    ggplot(aes(x=lambda_t, y=lambda_est, col=csp)) +
    geom_abline(slope=slope, intercept=0, col="black") +
    geom_point() +
    annotate("text",
      x=0, y=xymax,
      label=paste0("rho = ", cor_value),
      hjust=0.5, vjust=1, size=6) +
    xlab("Species loading targets") +
    ylab("Species loading estimates") +
    theme_bw(base_size=18) +
    theme(legend.position="none") +
    scale_color_manual(values=c(grey(0.6), "red")) +
    coord_fixed(ratio=1, xlim=c(-xymax, xymax),
                ylim=c(-xymax, xymax))
  return(p)
}

# End of file
