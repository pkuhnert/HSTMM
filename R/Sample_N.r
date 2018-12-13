Sample_N <- function(n, N, lambda, h_n, alpha_theta, beta_theta, j, m, nT){



  Nstar <- matrix(NA, nrow = m, ncol = nT+1)
  for(i in 1:m){
    for(k in 1:nT){
      Nstar[i,k+1] <- runifdisc(1, min=N[[j-1]][i, k+1] - h_n, max=N[[j-1]][i, k+1] + h_n)

    }
  }


  top1 <- top2 <- bot1 <- bot2 <- 0
  for(i in 1:m){
    for(k in 1:nT){
      # numerator

      top1 <- top1 + log(dpois(Nstar[i,k+1], lambda[[j]][i,k]))
      top2 <- dbetabin(k = n[i,k], n = Nstar[i,k+1], alpha = alpha_theta, beta = beta_theta)


      # denominator
      bot1 <- bot1 + log(dpois(N[[j-1]][i,k+1], lambda[[j]][i,k]))
      bot2 <- dbetabin(k = n[i,k], n = N[[j-1]][i,k+1], alpha = alpha_theta, beta = beta_theta)
    }
  }

  top <- top1 + top2
  bot <- bot1 + bot2


  pstar_N <- top - bot

  if(pstar_N > runif(1,0,1)){
    cat("New N\n")

    N[[j]] <- Nstar
    accept_N <- 1
  }
  else{
    N[[j]] <- N[[j-1]]
    accept_N <- 0
  }


  list(N = N, accept_N = accept_N)


}
