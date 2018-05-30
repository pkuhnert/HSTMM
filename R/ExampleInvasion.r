#' Example Invasion Script
#'
#'
#' @export



ExampleInvasion <- function(m = 6, nT = 10, niter = 100){

  # Setup
  longitude <- seq(40, 50, length = m)
  latitude <- seq(20, 30, length = m)

  coords <- data.frame(longitude = longitude, latitude = latitude)
  invasion_su <- SetupHSTMM(nT = 10, niter = 100, coords, m = m, tune = list(h_lambda = 0.1, h_tau = 0.3, h_K = 20, h_r = 0.1),
                         prior_lambda = list(lambda_m = log(15), lambda_sig = log(1.1)),
                         prior_tau = list(tau_m = log(2), tau_sig = log(1.1)),
                         prior_K = list(alpha = 40, beta = 50),
                         prior_r = list(mu = 0.5, sig2 = 0.25),
                         initialV = list(K = 150, r = 0.3))

  acceptR <- list(accept_lambda = NULL, accept_tau = NULL, accept_K = NULL, accept_r = NULL)


  # Run model
  MM_out <- HootenMM(niter = niter, m = m, nT = nT, r = invasion_su$r, lambda = invasion_su$lambda,
           K = invasion_su$K, coords = coords, tau = invasion_su$tau, N = invasion_su$N,
           prior_lambda = invasion_su$prior_lambda, prior_tau = invasion_su$prior_tau,
           prior_K = invasion_su$prior_K, prior_r = invasion_su$prior_r,
           acceptR = acceptR, Hparams = invasion_su$tune)


  MM_out

}
