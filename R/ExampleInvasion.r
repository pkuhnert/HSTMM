#' Example Invasion Script
#'
#' @param m number of spatial locations
#' @param nT number of time points
#' @param niter number of iterations
#' @param coords coordinates
#' @param prior_lambda prior for lambda
#' @param prior_tau prior for tau
#' @param prior_K prior for K
#' @param prior_r prior for r
#' @param prior_theta prior for theta
#' @param initialV initialV
#' @param tune tuning parameters
#' @param sample_n sample_n
#'
#'
#' @export



ExampleInvasion <- function(m = 6, nT = 10, niter = 100, coords,
                            prior_lambda = prior_lambda, prior_tau = prior_tau,
                            prior_K = prior_K, prior_r = prior_r, prior_theta,
                            initialV = initialV,
                            tune = list(h_lambda = 0.1, h_tau = 0.3, h_K = 20, h_r = 0.1,
                                        h_N = 1),
                            sample_n){

  # Setup

  invasion_su <- SetupHSTMM(nT = nT, niter = niter, coords, m = m,
                            tune = tune,
                         prior_lambda = prior_lambda,
                         prior_tau = prior_tau,
                         prior_K = prior_K,
                         prior_r = prior_r,
                         prior_theta = prior_theta,
                         initialV = initialV)

  acceptR <- list(accept_lambda = NULL, accept_tau = NULL, accept_K = NULL, accept_r = NULL)


  # Run model
  MM_out <- HootenMM(niter = niter, m = m, nT = nT, r = invasion_su$r, lambda = invasion_su$lambda,
           K = invasion_su$K, coords = coords, tau = invasion_su$tau, N = invasion_su$N,
           prior_lambda = invasion_su$prior_lambda, prior_tau = invasion_su$prior_tau,
           prior_K = invasion_su$prior_K, prior_r = invasion_su$prior_r,
           prior_theta = invasion_su$prior_theta,
           acceptR = acceptR, Hparams = invasion_su$tune, sample_n = sample_n)


  MM_out

}
