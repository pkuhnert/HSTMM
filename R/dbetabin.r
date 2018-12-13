dbetabin <- function(k, n, alpha, beta){
  
  term1 <- lgamma(n+1) + lgamma(k + alpha) + lgamma(n-k+beta) + 
    lgamma(alpha + beta)
  term2 <- lgamma(k+1) + lgamma(n-k+1) + lgamma(n+alpha+beta) +
    lgamma(alpha) + lgamma(beta)
  
  term1 - term2
  
  
}