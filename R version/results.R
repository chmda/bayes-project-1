source("R version/main_R.R")
source("R version/main_seq.R")

n <- 200

X <- MH_Gibbs_seq(n,init)

extract_theta <- function(X){
  n <- length(X)
  X_theta <-  c()
  for (i in 1:n){
    X_theta <- append(X_theta,X[[i]]$theta)
  }
return (X_theta)
}

extract_phi <- function(X){
  n <- length(X)
  X_phi <-  c()
  for (i in 1:n){
    X_phi <- append(X_phi,X[[i]]$phi)
  }
  return (X_phi)
}

extract_gamma_k <- function(X,k){
  n <- length(X)
  X_gamma_k <-  c()
  for (i in 1:n){
    X_gamma_k <- append(X_gamma_k,X[[i]]$gamma[k])
  }
  return (X_gamma_k)
}

extract_beta_k <- function(X,k){
  n <- length(X)
  X_beta_k <-  c()
  for (i in 1:n){
    X_beta_k <- append(X_beta_k,X[[i]]$beta[k])
  }
  return (X_beta_k)
}


results <- function(X){
  res <- list()
  res$beta_1 <- extract_beta_k(X,1)
  res$beta_2 <- extract_beta_k(X,2)
  res$beta_3 <- extract_beta_k(X,3)
  res$beta_4 <- extract_beta_k(X,4)
  res$beta_5 <- extract_beta_k(X,5)
  res$beta_6 <- extract_beta_k(X,6)
  res$beta_7 <- extract_beta_k(X,7)
  res$beta_8 <- extract_beta_k(X,8)

  res$gamma_1 <- extract_gamma_k(X,1)
  res$gamma_2 <- extract_gamma_k(X,2)
  res$gamma_3 <- extract_gamma_k(X,3)

  res$theta <- extract_theta(X)
  res$phi <- extract_phi(X)
  
  
  return (res)
}

res <- results(X)
mapply(summary,res)
