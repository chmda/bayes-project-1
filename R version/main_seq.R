source("R version/__init__.R")

library(rlang)
library(mvtnorm)
library(matlib)

# Hyper-paramètres 
prop_sd <- 0.01

Sigma <- inv(R) 
tau <- 0.001
mn = c(0,0,0)
prec <- diag(3)*100

#Initialisation

init <- list(theta = 0, phi = 0.01, gamma = c(0,0,0), beta = c(0,0,0,0,0,0,0,0),
             T = structure(.Data = c(10, -3, -3,
                                     -3, 135, -65,
                                     -3, -65, 135), .Dim = c(3, 3)),
             alpha = structure(.Data = c( 0,0,0,
                                          0,0,0,
                                          0,0,0,
                                          0,0,0,
                                          0,0,0,
                                          0,0,0,
                                          0,0,0,
                                          0,0,0,
                                          0,0,0,
                                          0,0,0,
                                          0,0,0,
                                          0,0,0,
                                          0,0,0,
                                          0,0,0,
                                          0,0,0,
                                          0,0,0,
                                          0,0,0,
                                          0,0,0,
                                          0,0,0,
                                          0,0,0,
                                          0,0,0,
                                          0,0,0,
                                          0,0,0,
                                          0,0,0,
                                          0,0,0,
                                          0,0,0,
                                          0,0,0,
                                          0,0,0,
                                          0,0,0,
                                          0,0,0,
                                          0,0,0,
                                          0,0,0,
                                          0,0,0,
                                          0,0,0,
                                          0,0,0,
                                          0,0,0,
                                          0,0,0,
                                          0,0,0), .Dim = c(38, 3))
) 

mu0 <- init$gamma 
T0 <- init$T

mu_tau <- function(i,param){ # Calcul de mu et tau pour un individu i
  k <- data$school[i]
  mu <- param$beta[1] * data$LRT[i]^2 + param$beta[2]*data$VR[i,2] + param$beta[3] * data$Gender[i] + param$beta[4] * data$School_gender[i,1] + param$beta[5] * data$School_gender[i,2] + param$beta[6] * data$School_denom[i,1] + param$beta[7] * data$School_denom[i,2] + param$beta[8] * data$School_denom[i,3] + param$alpha[k,1] + param$alpha[k,2] * data$LRT[i] + param$alpha[k,3] * data$VR[i,1]
  tau <- exp(param$theta + param$phi * data$LRT[i])
  return (c(mu,tau))
}

mu_tau_j <- function(j,param){ #Calcuk de lu et tau pour une école j
  mu <- param$beta[1] * data$LRT[data$school== j]^2 + param$beta[2]*data$VR[data$school== j,2] + param$beta[3] * data$Gender[data$school== j]
  + param$beta[4] * data$School_gender[data$school== j,1] + param$beta[5] * data$School_gender[data$school== j,2] 
  + param$beta[6] * data$School_denom[data$school== j,1] + param$beta[7] * data$School_denom[data$school== j,2]
  + param$beta[8] * data$School_denom[data$school== j,3]
  + param$alpha[j,1] + param$alpha[j,2] * data$LRT[data$school== j] + param$alpha[j,3] * data$VR[data$school== j,1]
  
  tau <- exp(param$theta + param$phi * data$LRT[data$school== j])
  
  return (c(mu,tau))
}


MH_Gibbs_seq <- function(n,init){
  X <- rep(list(init),n+1)
  X[[1]] <- init
  
  for (l in 2:(n+1)){
    X[[l]] <- duplicate(X[[l-1]])

    # Calcul des beta 
    for (k in 1:8){
      nu_k <- switch(k, # Variable associée à beta_k
                     data$LRT^2,data$VR[,2] ,data$Gender
                     ,  data$School_gender[,1] , data$School_gender[,2] 
                     ,  data$School_denom[,1] ,  data$School_denom[,2]
                     ,  data$School_denom[,3])
      
      tau_i <- 1:N 
      mu_i <- 1:N
      a_i <- 1:N
      for (i in 1:N){  
        L <- mu_tau(i,X[[l]])
        mu_i <- L[1]
        tau_i[i] <- L[2]
        a_i[i] <- data$Y[i] - mu_i + nu_k[i]*X[[l]]$beta[k]
      }
      sigma2 <- 1/(tau + sum(nu_k^2*tau_i))
      mu_2 <- sigma2*(sum(nu_k*tau_i*a_i))  
      X[[l]]$beta[k] <- rnorm(1,mu_2,sqrt(sigma2))
    }
    
    # Calcul des alpha 
    
    for (j in 1:M){ # On échantillone alpha_j pour chaque école j
      prop <- X[[l]]$alpha[j,] + rnorm(3, 0, sd=prop_sd)
      mat_T <- X[[l]]$T
      gamma <- X[[l]]$gamma

      
      
      A_2 <- data$LRT[data$school == j]
      A_1 <- rep(1, length(A_2))
      A_3 <- data$VR[data$school == j,1]
      Aj <- matrix(NA,3,length(A_2))
      Aj[1,] <- A_1
      Aj[2,] <- A_2
      Aj[3,] <- A_3
      
      L <- mu_tau_j(j,X[[l]])
      mu_j <- L[1]
      tau_j <- L[2]
      
      log_pdf <- function(alpha){
        p1 <- -0.5*mahalanobis(alpha,gamma,mat_T,inverted=TRUE)
        
        a_j <- data$Y[data$school == j] - mu_j + alpha%*%Aj
        
        p2 <- -0.5 * sum( ( (alpha%*%Aj)^2 - 2*a_j*(alpha%*%Aj) )*tau_j)
        
        return (p1 + p2)
      }
      
      top <- log_pdf(prop)
      bottom <- log_pdf(X[[l]]$alpha[j,])
      acc <- exp(top - bottom)
    
      
      if (runif(1) < acc) {
        X[[l]]$alpha[j,] <- prop 
      }
    }
    
    # theta
    prop <- X[[l]]$theta + rnorm(1, 0, sd=prop_sd)
    tau_i <- 1:N
    mu_i <- 1:N
    for (i in 1:N){
      L <- mu_tau(i,X[[l]])
      mu_i[i] <- L[1]
      tau_i[i] <- L[2]
    }
    
    phi <- X[[l]]$phi
    log_pdf <- function(theta){
      -0.5*theta^2*tau - 0.5 * (data$N*theta + sum( (data$Y - mu_i)^2*exp(theta + phi * data$LRT)))
    }
    
    top <- log_pdf(prop)
    bottom <- log_pdf(X[[l]]$theta)
    acc <- exp(top - bottom)
    
    if (runif(1) < acc) {
      X[[l]]$theta <- prop 
    }
    
    # phi
    prop <- X[[l]]$phi + rnorm(1, 0, sd=prop_sd)
    mu_i <- 1:N
    for (i in 1:N){
      L <- mu_tau(i,X[[l]])
      mu_i[i] <- L[1]
    }
    
    theta <- X[[l]]$theta
    log_pdf <- function(phi){
      -0.5*phi^2*tau - 0.5 * sum(phi*data$LRT +  (data$Y - mu_i)^2*exp(theta + phi*data$LRT))
    }
    
    top <- log_pdf(prop)
    bottom <- log_pdf(X[[l]]$phi)
    acc <- exp(top - bottom)
    
    
    if (runif(1) < acc) {
      X[[l]]$phi <- prop 
    }
    
    # gamma
    prec <- T0 + data$M*X[[l]]$T
    alpha_sum <- c(sum(X[[l]]$alpha[,1]),sum(X[[l]]$alpha[,2]),sum(X[[l]]$alpha[,3]))
    mu <- solve(a=prec,b=(T0 %*% mu0 + X[[l]]$T %*% alpha_sum ))
    X[[l]]$gamma <- as.vector(rmvnorm(1,mu, inv(prec)))
      
    # Calcul de T
    
    S <- matrix(0,3,3)
    for (j in 1:data$M){
 
      S <- S + ( (X[[l]]$alpha[j,] - X[[l]]$gamma) %*% t(X[[l]]$alpha[j,] - X[[l]]$gamma) )
    }
    R_new <- inv( inv(R) + S)
    X[[l]]$T <- array(rWishart(1,data$M + 3,R_new), dim = c(3, 3)) # 3x3x1 donc on redimensionne
    
    
  }
  
  return (X)
}