#' Sample the dispersal parameter, tau
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
#' @param h_tau M-H step value for tau
#' @param h_lambda M-H step value for lambda
#' @param mu_tau prior mean for tau
#' @param sigma_tau prior variance matrix for tau
#'
#' @description Simulates from the posterior distribution of the dispersal parameter, tau.
#'
#' @import mvtnorm
#' @import stats
#'
#' @export




Sample_tau <- function(tau, j, m, nT, r, lambda, K, coords, N, h_tau, h_lambda, mu_tau, sigma_tau){

  longitude <- coords$longitude
  latitude <- coords$latitude

  Ltaustar <- rnorm(m, log(tau[[j-1]]), h_tau)
  taustar <- exp(Ltaustar)


 lambdastar <- Llambdastar  <- matrix(NA, nrow = m, ncol = nT)
 Llambdastar[,1] <- rnorm(m, log(lambda[[j]][,1]), h_lambda)
 lambdastar[,1] <- exp(Llambdastar[,1])


 G <- vector("list", length = nT)
 G[[1]] <- Create_G(m = m, r = r[j-1], lambda = lambda[[j]], K = K[j-1], t = 1)
 M <- Create_M(m = m, longitude = longitude, latitude = latitude, tau = taustar)

 # Run forward
 for(t in 2:(nT)){
   lambdastar[,t] <- M %*% G[[t-1]] %*% lambdastar[,t-1]
   G[[t]] <- Create_G(m = m, r = r[j-1], lambda = lambdastar, K = K[j-1], t = t)

 }


   term1 <- term2 <- NULL

   # proposed value
   for(t in 1:nT){
      term1[t] <- sum(log(dpois(N[[j]][,t+1], lambdastar[,t])))
      term2[t] <- sum(log(dpois(N[[j]][,t+1], lambda[[j]][,t])))
   }


   top <- sum(term1) + sum(dmvnorm(log(taustar),mu_tau, sigma_tau, log = TRUE))
   bot <- sum(term2) + sum(dmvnorm(log(tau[[j-1]]), mu_tau, sigma_tau, log = TRUE))

 pstar_tau <- top-bot

 if(pstar_tau > runif(1,0,1)){
   tau[[j]] <- taustar
   accept_tau <- 1
 }
 else{
   tau[[j]] <- tau[[j-1]]
   accept_tau <- 0
 }


 list(tau = tau, accept_tau = accept_tau)



}
