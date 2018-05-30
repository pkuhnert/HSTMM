#' Sample the Carrying Capacity, K
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
#' @param h_K M-H step value for K
#' @param h_lambda M-H step value for lambda
#' @param alpha_K prior parameter for K
#' @param beta_K prior parameter for K
#'
#' @description Simulates from the posterior distribution of the carrying capacity parameter, K.
#'
#' @export

Sample_K <- function(tau, j, m, nT, r, lambda, K, coords, N, h_K, h_lambda, alpha_K, beta_K){

  Kstar <- rnorm(1, K[j-1], h_K)
#  cat("K = ", K[j-1], "  Kstar = ", Kstar, "\n")

  lambdastar <- Llambdastar  <- matrix(NA, nrow = m, ncol = nT)
  Llambdastar[,1] <- rnorm(m, log(lambda[[j]][,1]), h_lambda)
  lambdastar[,1] <- exp(Llambdastar[,1])

  G <- vector("list", length = nT)
  G[[1]] <- Create_G(m = m, r = r[j-1], lambda = lambda[[j]], K = Kstar, t = 1)
  M <- Create_M(m = m, longitude = longitude, latitude = latitude, tau = tau[[j]])

  # Run forward
  for(t in 2:(nT)){
    lambdastar[,t] <- M %*% G[[t-1]] %*% lambdastar[,t-1]
    G[[t]] <- Create_G(m = m, r = r[j-1], lambda = lambdastar, K = Kstar, t = t)

  }


  term1 <- term2 <- NULL

  # proposed value
  for(t in 1:nT){
    term1[t] <- sum(log(dpois(N[[j]][,t+1], lambdastar[,t])))
    term2[t] <- sum(log(dpois(N[[j]][,t+1], lambda[[j]][,t])))
  }


  top <- sum(term1) + sum(dgamma(Kstar, alpha_K, beta_K))
  bot <- sum(term2) + sum(dgamma(K[j-1], alpha_K, beta_K))

  pstar_K <- top - bot

  if(pstar_K > runif(1,0,1)){
    cat("New K\n")
    K[j] <- Kstar
    accept_K <- 1
  }
  else{
    K[j] <- K[j-1]
    accept_K <- 0
  }


  list(K = K, accept_K = accept_K)



}
