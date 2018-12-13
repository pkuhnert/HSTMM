#' Sample Growth Rate, r
#'
#' @param tau dispersal parameter (vector of length m)
#' @param j MC iteration step
#' @param m number of spatial locations
#' @param nT number of time points
#' @param r growth rate (scalar)
#' @param lambda intensity matrix (m x nT)
#' @param K carrying capacity (scalar)
#' @param coords matrix of spatial co-ordinates of locations (longitude, latitude)
#' @param N true number of organisms (m x nT)
#' @param h_r M-H step value for r
#' @param h_lambda M-H step value for lambda
#' @param mu_r prior mean for r
#' @param sig2_r prior variance for r
#'
#' @description Simulates from the posterior distribution of the growth rate, r.
#'
#' @import stats
#'
#' @export


Sample_r <- function(tau, j, m, nT, r, lambda, K, coords, N, h_r, h_lambda, mu_r, sig2_r){

  longitude <- coords$longitude
  latitude <- coords$latitude

  # proposal
  rstar <- rnorm(1, r[j-1], h_r)

  lambdastar <- Llambdastar  <- matrix(NA, nrow = m, ncol = nT)
  Llambdastar[,1] <- rnorm(m, log(lambda[[j]][,1]), h_lambda)
  lambdastar[,1] <- exp(Llambdastar[,1])

  G <- vector("list", length = nT)
  G[[1]] <- Create_G(m = m, r = rstar, lambda = lambda[[j]], K = K[j], t = 1)
  M <- Create_M(m = m, longitude = longitude, latitude = latitude, tau = tau[[j]])

  # Run forward
  for(t in 2:(nT)){
    lambdastar[,t] <- M %*% G[[t-1]] %*% lambdastar[,t-1]
    G[[t]] <- Create_G(m = m, r = rstar, lambda = lambdastar, K = K[j], t = t)

  }


  term1 <- term2 <- NULL

  # proposed value
  for(t in 1:nT){
    term1[t] <- sum(log(dpois(N[[j]][,t+1], lambdastar[,t])))
    term2[t] <- sum(log(dpois(N[[j]][,t+1], lambda[[j]][,t])))
  }


  top <- sum(term1) + sum(dnorm(rstar, sig2_r))
  bot <- sum(term2) + sum(dnorm(r[j-1], mu_r, sqrt(sig2_r)))



  pstar_r <- top - bot

  if(pstar_r > runif(1,0,1)){
    r[j] <- rstar
    accept_r <- 1
  }
  else{
    r[j] <- r[j-1]
    accept_r <- 0
  }


  list(r = r, accept_r = accept_r)



}
