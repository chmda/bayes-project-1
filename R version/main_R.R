source("R version/__init__.R")

library(mvtnorm)

mu_tau <- function(i,param){
  mu <- param[1] * data$LRT[i]^2 + param[2]*data$VR[i][2] + param[3] * data$Gender[i]
         + param[4] * data$School_gender[i][1] + param[5] * data$School_gender[i][2] 
         + param[6] * data$School_denom[i][1] + param[7] * data$School_denom[i][2]
         + param[8] * data$School_denom[i][3]
         + param[9] + param[10] * data$LRT[i] + param[11] * data$VR[i][1]
  
  tau <- exp(param[12] + param[13] * data$LRT[i])
  
  return (c(mu,tau))
}

mu_tau_j <- function(j,param){
  mu <- param[1] * data$LRT[data$school== j]^2 + param[2]*data$VR[data$school== j][2] + param[3] * data$Gender[data$school== j]
  + param[4] * data$School_gender[data$school== j][1] + param[5] * data$School_gender[data$school== j][2] 
  + param[6] * data$School_denom[data$school== j][1] + param[7] * data$School_denom[data$school== j][2]
  + param[8] * data$School_denom[data$school== j][3]
  + param[9] + param[10] * data$LRT[data$school== j] + param[11] * data$VR[data$school== j][1]
  
  tau <- exp(param[12] + param[13] * data$LRT[data$school== j])
  
  return (c(mu,tau))
}

prop_sd <- 0.01
Sigma <- inv(R) 
mn = c(0,0,0)
prec <- diag(3)*100

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



sample_beta <- function(k,X){
  nu_k <- switch(k, # Variable associée à beta_k
                 data$LRT^2,data$VR[,2] ,data$Gender
                 ,  data$School_gender[,1] , data$School_gender[,2] 
                 ,  data$School_denom[,1] ,  data$School_denom[,2]
                 ,  data$School_denom[,3])
  tau_i <- 1:N
  a_i <- 1:N
  for (i in 1:N){
    L <- mu_tau(i,X)
    mu <- L[1]
    tau_i[i] <- L[2]
    nu_ki <- nu_k[i]
    a_i[i] <- data$Y[i] - mu + nu_ki*X[k]
  }
  va <- 1/(tau + sum(nu_k*tau_i))
  mu_2 <- va*(sum(nu_k*tau_i - a_i)) 
  X[k] <- rnorm(mu_2,sqrt(va))
  return (X)
}


sample_alpha <-function(X){ # Marche aléatoire 
  alpha <- matrix(NA,M,3)
  for (j in 1:M){
    prop <- X[9:11][,j] + rnorm(3, 0, sd=prop_sd)
    mat_T <- rWishart(1,3,Sigma)
    gamma <- rmvnorm(1, mn, prec)
    A_2 <- data$LRT[data$school == j]
    A_1 <- rep(1, length(A_2))
    A_3 <- data$VR[data$school == j][,1]
    Aj <- c(A_1,A_2,A_3)
    L <- mu_tau_j(j,X)
    mu_j <- L[1]
    tau_j <- L[2]
    log_pdf <- function(alpha){
      p1 <- -0.5*t(alpha-gamma)%*%T%*%(alpha-gamma)
      b_j <- mu_j - alpha%*%Aj - data$Y[data$school == j]
      p2 <- -0.5 * sum( (alpha%*%Aj)^2 - 2*b_j%*%(alpha%*%Aj)*tau_j)
      return (p1 + p2)
    }
    
    top <- log_pdf(prop)
    bottom <- log_pdf(X[9:11][,j])
    acc <- exp(top - bottom)
    
    if (runif() < acc) {
      alpha[j] <- prop 
    }
    else {
      alpha[j] <- X[9:11][,j]
    }
  }
  X[9:11] <- t(alpha)
  return (X)
}

sample_theta <- function(X){
  prop <- X[12] + rnorm(1, 0, sd=prop_sd)
  tau_i <- 1:N
  mu_i <- 1:N
  for (i in 1:N){
    L <- mu_tau(i,X)
    mu_i[i] <- L[1]
    tau_i[i] <- L[2]
  }
  
  log_pdf <- function(theta){
    -0.5*theta^2*tau - 0.5 * (theta + sum( (data$Y - mu_i)^2*exp(theta)*tau_i))
  }
  
  top <- log_pdf(prop)
  bottom <- log_pdf(X[12])
  acc <- exp(top - bottom)
  
  if (runif() < acc) {
    X[12] <- prop 
  }
  
  return (X)
}

sample_phi <- function(X){
  prop <- X[13] + rnorm(1, 0, sd=prop_sd)
  tau_i <- 1:N
  mu_i <- 1:N
  for (i in 1:N){
    L <- mu_tau(i,X)
    mu_i[i] <- L[1]
    tau_i[i] <- L[2]
  }
  
  log_pdf <- function(phi){
    -0.5*phi^2*tau - 0.5 * sum(phi*data$LRT +  (data$Y - mu_i)^2*exp(phi)*data$LRT*tau_i)
  }
  
  top <- log_pdf(prop)
  bottom <- log_pdf(X[13])
  acc <- exp(top - bottom)
  
  if (runif() < acc) {
    X[13] <- prop 
  }
  
  return (X)
}

MH_Gibbs <- function(n,init){
  X <- 1:n
  X[1] <- init
  for (i in 1:n){
    for (k in 1:8){
      X[i] <- sample_beta(k,X)
    }
    X[i] <- sample_alpha(X)
    X[i] <- sample_theta(X)
    X[i] <- sample_phi(X)
  }
  return (X)
}
