require(coda)
require(ggplot2)

##' @title Transform an array into an MCMC
##' @param x An array object
##' @return An mcmc object
##' @author Ghislain Vieilledent
arr2mcmc <- function(x) {
  return(mcmc(as.data.frame(x),
              start=burnin+1 , end=ngibbs, thin=thin))
}

##' @title 
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
  nsamp <- nrow(mcmc_list_alpha[[1]])

  # psrf gelman indice 
  psrf_alpha <- max(gelman.diag(mcmc_list_alpha,
                                 multivariate=FALSE)$psrf[,2])
  psrf_V_alpha <- gelman.diag(mcmc_list_V_alpha)$psrf[,2]
  psrf_beta <- max(gelman.diag(mcmc_list_sp[, grep("beta", colnames(mcmc_list_sp[[1]]))],
                                multivariate=FALSE)$psrf[,2])
  psrf_lambda <- max(gelman.diag(mcmc_list_lambda,
                                  multivariate=FALSE)$psrf[,2], na.rm=TRUE)
  psrf_lv <- max(gelman.diag(mcmc_list_lv,
                              multivariate=FALSE)$psrf[,2 ])

  Rhat <- data.frame(Variable=c("alpha", "Valpha", "beta", "lambda", "W"),
                     maxRhat=c(psrf_alpha, psrf_V_alpha, psrf_beta, psrf_lambda, psrf_lv))

  return(Rhat)

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
##'   "lambdas" for all lambdas. Could also be
##'   "(sp_1\\.lambda_1|sp_4\\.lambda_1)" for example.
##' @return
##' @author Ghislain Vieilledent
mcmc_lambdas <- function(model_output, re="lambdas") {
  # Parameters for coda object
  burnin <- model_output[[1]]$model_spec$burnin
  ngibbs <- burnin + model_output[[1]]$model_spec$mcmc
  thin <-  model_output[[1]]$model_spec$thin
  # MCMC lists
  mcmc_list_sp <- mcmc.list(lapply(lapply(model_output,"[[","mcmc.sp"), arr2mcmc))
  mcmc_list_lambda <- mcmc.list(
    lapply(mcmc_list_sp[, grep(re, colnames(mcmc_list_sp[[1]]), value=TRUE)], arr2mcmc))
  return(mcmc_list_lambda)
}

# End of file
