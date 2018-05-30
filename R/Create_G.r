#' Generate Growth Matrix
#'
#' @param m number of spatial locations
#' @param r growth rate parameter
#' @param lambda intensity matrix (m x nT)
#' @param K carrying capacity (scalar)
#' @param t time point
#' @param plot (Default: FALSE)
#'
#' @description Computes the growth matrix, G. This represents a diagonal
#' matrix indexed by t.

Create_G <- function(m, r, lambda, K, t, plot = FALSE){


    G <- diag(m)
    for(i in 1:m)
      G[i,i] <- exp(r * (1-lambda[,t][i]/K))



  if(plot){
      image(G, col = rev(heat.colors(12)), zlim = c(0, 2), main = paste("Time: ", t))
      box()
  }

  G
}
