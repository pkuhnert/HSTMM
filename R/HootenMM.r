#' Hooten Space-Time Matrix Model
#'
#' @param niter number of iterations
#' @param m number of spatial locations
#' @param nT number of time points
#' @param r growth rate (scalar)
#' @param lambda intensity matrix (m x nT)
#' @param K carrying capacity (scalar)
#' @param coords matrix of spatial co-ordinates of locations (longitude, latitude)
#' @param tau dispersal parameter (vector of length m)
#' @param N true number of organisms (m x nT)
#' @param acceptR acceptance rates of parameters for M-H step
#' @param Hparams hyper-parameters for prior probability distributions
#' @param prior_lambda prior for lambda (mu_lambda, Sigma_lambda)
#' @param prior_tau prior for dispersal parameter, tau
#' @param prior_K prior for carrying capacity, K
#' @param prior_r prior for growth rate, r
#' @param prior_theta prior for theta
#' @param sample_n sample_n
#'
#' @description Uses the M-H algorithm to implement a BHM S-T model for invasions
#'
#' @details This implements the M-H algorithm outlined in Hooten et al. (2007).
#'
#' @references Hooten, M.B., Wikle, C.K., Dorazio, R.M. and Royle, J.A. (2007) Hierarchical Spatiotemporal Matrix
#' Models for Characterizing Invasions, Biometrics, 63, 558-567.
#'
#' @import stats
#' @importFrom svMisc progress
#'
#' @export


HootenMM <- function(niter, m, nT, r, lambda, K, coords, tau, N, acceptR, Hparams,
                     prior_lambda, prior_tau, prior_K, prior_r, prior_theta, sample_n){


  #----------------------- Extraction ------------------------#

  # spatial co-ordinates
  longitude <- coords$longitude
  latitude <- coords$latitude

  # acceptance rates
  accept_lambda <- acceptR$lambda
  accept_tau <- acceptR$tau
  accept_K <- acceptR$K
  accept_r <- acceptR$r

  # Hyper-parameters
  h_lambda <- Hparams$h_lambda
  h_tau <- Hparams$h_tau
  h_K <- Hparams$h_K
  h_r <- Hparams$h_r
  h_n <- Hparams$h_n



  accept_lambda <- accept_tau <- accept_K <- accept_r <- accept_N <- NULL


  #----------------------- Initialisation ------------------------#
  # initial values for the process at j=1
  j <- 1
  progress(j, progress.bar = TRUE)
  G <- vector("list", length = nT)
  G[[1]] <- Create_G(m = m, r = r[j], lambda = lambda[[j]], K = K[j], t = 1, plot = FALSE)
  M <- Create_M(m = m, longitude = longitude, latitude = latitude, tau = tau[[j]], plot = FALSE)


  # run forward process model to obtain initial values of the process (#4 in Biometrics paper)
  for(t in 2:(nT+1)){
    lambda[[j]][,t] <- M %*% G[[t-1]] %*% lambda[[j]][,t-1]

    N[[j]][,t] <- rpois(m, lambda[[j]][,t])
    G[[t]] <- Create_G(m = m, r = r[j], lambda = lambda[[j]], K = K[j], t = t, plot = FALSE)

  }


  #----------------------- Sampling ------------------------#

  cat("MCMC Run \n")
  for(j in 2:niter){
    progress(j, progress.bar = TRUE)

    # Sample lambda
    val_lambda <- Sample_lambda(m = m, nT = nT, r = r, j = j, lambda = lambda, K = K, coords = coords, tau = tau, N = N,
                                mu_lambda = prior_lambda$mu_lambda, sigma_lambda = prior_lambda$Sigma_lambda, h_lambda = h_lambda)
   # N <- val_lambda$N
    lambda <- val_lambda$lambda
    accept_lambda[j-1] <- val_lambda$accept_lambda

    # Sample N
    val_N <- Sample_N(n = sample_n, N = N, lambda = lambda, h_n = h_n,
                      alpha_theta = prior_theta$alpha_theta,
                      beta_theta = prior_theta$beta_theta, j = j, m = m, nT = nT)
    N <- val_N$N
    accept_N[j-1] <- val_N$accept_N



    # sample tau
    val_tau <- Sample_tau(tau = tau, j = j, m = m, nT = nT, r = r, lambda = lambda, K = K, coords = coords,
                          N = N, h_tau = h_tau, h_lambda = h_lambda, mu_tau = prior_tau$mu_tau, sigma_tau = prior_tau$Sigma_tau)

    tau <- val_tau$tau
    accept_tau[j-1] <- val_tau$accept_tau


    # Sample K
    val_K <- Sample_K(tau = tau, j = j, m = m, nT = nT, r = r, lambda = lambda, K = K, coords = coords, N = N, h_K = h_K,
                      h_lambda = h_lambda, alpha_K = prior_K$alpha, beta_K = prior_K$beta)
    K <- val_K$K
    accept_K[j-1] <- val_K$accept_K

    # Sample r
    val_r <- Sample_r(tau = tau, j = j, m = m, nT = nT, r = r, lambda = lambda, K = K, coords = coords, N = N, h_r = h_r,
                      h_lambda = h_lambda, mu_r = prior_r$mu, sig2_r = prior_r$sig2)
    r <- val_r$r
    accept_r[j-1] <- val_r$accept_r

    if(j == niter) cat("Done!\n")

  }

  AR <- list(accept_lambda = accept_lambda, accept_tau = accept_tau,
             accept_K = accept_K, accept_r = accept_r, accept_N = accept_N)


  list(N = N, lambda = lambda, tau = tau, K = K, r = r, AR = AR)


}



