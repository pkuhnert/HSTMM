#' Generate Movement Matrix
#'
#' @param m number of spatial locations
#' @param longitude longitude of location
#' @param latitude latitude of location
#' @param tau dispersal parameter (vector of length m)
#' @param plot (Default: FALSE)
#'
#' @importFrom stats "dist"
#' @import graphics
#' @import grDevices
#'
#' @description Computes the movement matrix, M, based on a spatially varying Gaussian kernel

Create_M <- function(m, longitude, latitude, tau, plot = FALSE){

  M <- matrix(0, nrow = m, ncol = m)
  for(i in 1:m){
    for(k in 1:m){
      M[i,k] <- exp(-(dist(matrix(c(longitude[i], latitude[i], longitude[k], latitude[k]),
                                  byrow = TRUE, ncol = 2))^2)/tau[k])
    }
  }


  M <- apply(M, 1, function(x) x/sum(x))

  if(plot){
      image(M, col = rev(heat.colors(12)), zlim = c(0, 2) )
      box()

  }

  M


}


