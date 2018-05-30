#' Set up the HSTMM
#'
#' @param nT number of time steps
#' @param niter number of MCMC iterations
#' @param coords matrix of longitude and latitude co-ordinates
#' @param m number of spatial locations
#' @param tune tuning parameters for MH step. A value is listed for each type of parameter.
#' @param prior_lambda prior for lambda
#' @param prior_tau prior for tau
#' @param prior_K prior for K
#' @param prior_r prior for r
#' @param initialV initial values
#'
#' @import mvtnorm
#' @import stats
#'
#' @export
#'
#'




SetupHSTMM <- function(nT = 10, niter = 100, coords, m = nrow(coords),
                       tune = list(h_lambda = 0.1, h_tau = 0.3, h_K = 20, h_r = 0.1),
                       prior_lambda = list(lambda_m = log(15), lambda_sig = log(1.1)),
                       prior_tau = list(tau_m = log(2), tau_sig = log(1.1)),
                       prior_K = list(alpha = 40, beta = 50),
                       prior_r = list(mu = 0.5, sig2 = 0.25),
                       initialV = list(K = 150, r = 0.3)){


  longitude <- coords$longitude
  latitude <- coords$latitude

  # hyper-parameters
  mu_lambda <- rep(prior_lambda$lambda_m, m)
  mu_tau <- rep(prior_tau$tau_m, m)

  Sigma_lambda <- Sigma_tau <- matrix(NA, nrow = m, ncol = m)
  for(i in 1:m){
    for(k in 1:m){
      Sigma_lambda[i,k] <- prior_lambda$lambda_sig * exp(-2 * dist(matrix(c(longitude[i], latitude[i], longitude[k], latitude[k]),
                                                           byrow = TRUE, ncol = 2)))
      Sigma_tau[i,k] <- prior_tau$tau_sig * exp(-2 * dist(matrix(c(longitude[i], latitude[i], longitude[k], latitude[k]),
                                                               byrow = TRUE, ncol = 2)))

    }
  }


  # Set up parameters

  K <- r <- vector("numeric", length = niter)
  tau <- lambda <- N <- vector("list", length = niter)
  K[1] <- initialV$K
  r[1] <- initialV$r
  lambda[[1]] <- N[[1]] <- matrix(NA, nrow = m, ncol = nT+1)
  lambda[[1]][,1] <- rmvnorm(1, mu_lambda, sigma = Sigma_lambda)

  tau[[1]]<- rmvnorm(1, mu_tau, sigma = Sigma_tau)


  list(tune = tune, K = K, r = r, lambda = lambda, tau = tau, N = N,
       prior_lambda = list(mu_lambda = mu_lambda, Sigma_lambda = Sigma_lambda),
       prior_tau = list(mu_tau = mu_tau, Sigma_tau = Sigma_tau), prior_K = prior_K,
       prior_r = prior_r)



}

