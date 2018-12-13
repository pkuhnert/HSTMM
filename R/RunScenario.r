#' RunScenario
#'
#' @description Just a simple scenario to run Hooten Model
#'
#' @param niter Number of iterations
#' @param nT Number of time points
#' @param nburn Number of burn-ins
#'
#' @import ggplot2
#' @import mvtnorm
#' @import RColorBrewer
#' @import gridExtra
#'
#'
#' @export


RunScenario <- function(niter = 30000, nburn = 10000, nT = 10){


  if(nburn > niter) stop("No. of burn-ins needs to be less than number of iterations.\n")

  # create dummy locations
  longitude <- seq(40, 50, length = 5)
  latitude <- seq(20, 30, length = 5)
  coords <- expand.grid(longitude = longitude, latitude = latitude)
  m <- nrow(coords)

  coords_farm <- coords
  coords_farm$farm <- as.factor(c(rep(1, 10), rep(2, 15)))
  p_farm <- ggplot(coords_farm, aes_string("longitude", "latitude", col = "farm")) +
    geom_point(shape = 15, size = 5) +
    theme(axis.text=element_blank(),
          axis.ticks=element_blank(),
          axis.title=element_blank()) + coord_fixed() +
    geom_hline(yintercept = 23.75, size = 1, color = "black")

  lambda_m <- c(rep(log(15), 10), rep(log(1.5), 15))

  # Surveillance
  sample_n <- matrix(0, nrow = m, ncol = nT)


  Hout <- ExampleInvasion(m = m, nT = nT, niter = niter, coords = coords,
                         prior_lambda = list(lambda_m = lambda_m, lambda_sig = log(1.01)),
                         prior_tau = list(tau_m = log(1.5), tau_sig = log(1.01)),
                         prior_K = list(alpha = 3, beta = 12),
                         prior_r = list(mu = 0.5, sig2 = 0.25),
                         prior_theta = list(alpha_theta = 3, beta_theta = 12),
                         initialV = list(K = 150, r = 0.3),
                         tune = list(h_lambda = 0.05, h_tau = 0.2, h_K = 10, h_r = 0.05,
                                     h_n = 1),
                         sample_n = sample_n)

  Hres <- lapply(Hout$AR, mean)


  # iteration range to examine
  startMC <- niter-nburn+1
  stopMC <- niter


  #------------------- lambda ---------------------------#

  sumlambda <- Hout$lambda[[startMC]]
  count <- 0
  for(i in (startMC+1):stopMC){
    check <- (i/2) - trunc(i/2)
    if(check == 0){
      sumlambda <- sumlambda + Hout$lambda[[i]]
      count <- count + 1
    }

  }
  #lambda_m <- sumlambda/(stopMC - startMC + 1)
  lambda_m <- sumlambda/count
  lambda_m <- data.frame(lambda_m)
  names(lambda_m) <- paste("t", 1:nT, sep = "")

  coords_lambda <- cbind(coords, lambda_m)

  p <- list()
  for(i in 1:nT){
    sub <- coords_lambda[,c("longitude", "latitude", paste("t", i, sep = ""))]
    names(sub)[3] <- "time"
    p[[i]] <- ggplot(sub, aes(longitude, latitude, col = time)) +
      geom_point(size = 10, shape = 15) + ggtitle(paste("Time: ", i, sep = "")) +
      scale_color_gradientn("lambda (Intensity)",
                            colours=c(rev(brewer.pal(8,"Spectral"))), na.value = "white",
                            limits = range(0, lambda_m, 80)) +
      theme(axis.text=element_blank(),
            axis.ticks=element_blank(),
            axis.title=element_blank()) + coord_fixed()
  }



  p_spread <- marrangeGrob(p, nrow = 2, ncol = 2, as.table = FALSE)



 list(p_farm = p_farm, p_spread = p_spread, Hres = Hres)

}




