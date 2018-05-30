#' Sample the intensity, lambda
#'
#' @param m number of spatial locations
#' @param nT number of time points
#' @param r growth rate (scalar)
#' @param j MC iteration step
#' @param lambda intensity matrix (m x nT)
#' @param K carrying capacity (scalar)
#' @param coords matrix of spatial co-ordinates of locations (longitude, latitude)
#' @param tau dispersal parameter (vector of length m)
#' @param N true number of organisms (m x nT)
#' @param mu_lambda prior mean vector for lambda
#' @param sigma_lambda prior variance matrix for lambda
#' @param h_lambda M-H step value for lambda
#'
#' @description Simulates from the posterior distribution of the intensity, lambda.
#'
#' @import mvtnorm
#' @import stats
#'
#' @export

Sample_lambda <- function(m, nT, r, j, lambda, K, coords, tau, N,
                          mu_lambda, sigma_lambda, h_lambda){

  longitude <- coords$longitude
  latitude <- coords$latitude


  # sample from proposal
  lambdastar <- Llambdastar  <- matrix(NA, nrow = m, ncol = nT)
  Nstar <- matrix(NA, nrow = m, ncol = nT)
  Llambdastar[,1] <- rnorm(m, log(lambda[[j-1]][,1]), h_lambda)
  lambdastar[,1] <- exp(Llambdastar[,1])
  Nstar[,1] <- rpois(m, lambdastar[,1])

  G <- vector("list", length = nT)
  G[[1]] <- Create_G(m = m, r = r[j-1], lambda = lambdastar, K = K[j-1], t = 1)
  M <- Create_M(m = m, longitude = longitude, latitude = latitude, tau = tau[[j-1]])

  # Run forward
  for(t in 2:(nT)){
    lambdastar[,t] <- M %*% G[[t-1]] %*% lambdastar[,t-1]
    Nstar[,t] <- rpois(m, lambdastar[,t])
    G[[t]] <- Create_G(m = m, r = r[j-1], lambda = lambdastar, K = K[j-1], t = t)

  }


  term1 <- term3 <-  NULL
  for(t in 1:nT){
    term1[t] <- sum(log(dpois(N[[j-1]][,t+1], lambdastar[,t])))
    term3[t] <- sum(log(dpois(N[[j-1]][,t+1], lambda[[j-1]][,t])))
  }


  top <- sum(term1) + dmvnorm(log(lambdastar[,1]), mu_lambda, sigma_lambda, log = TRUE)
  bot <- sum(term3) + dmvnorm(log(lambda[[j-1]][,1]), mu_lambda, sigma_lambda, log = TRUE)
  pstar_lambda <- top - bot


  if(pstar_lambda > runif(1,0,1)){
    cat("New lambda\n")
    lambda[[j]] <- lambdastar
    N[[j]] <- cbind(rep(NA, m), Nstar)
    accept_lambda <- 1
  }
  else{
    lambda[[j]] <- lambda[[j-1]]
    N[[j]] <- N[[j-1]]
    accept_lambda <- 0
  }


  list(lambda = lambda, N = N, accept_lambda = accept_lambda)



}
