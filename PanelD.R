rm(list=ls()) ## clean workspace
# load necessary packages
library(tidyverse)
library(mvtnorm)
library(mgcv)
library(magrittr)
library(coneproj)  
library(gtools)
library(dplyr)
library(geometry)
library(uniformly)
library(pracma)
library(deSolve)
######## FUNCTIONS ##################


# function that computes the normalized feasibility from an interaction matrix analytically
# inputs: alpha = interaction matrix
#         individual = TRUE: community-level Omega
#                      FALSE: individual-level omega
# output: out = the normalized feasibility
# !!! double check community or species level
Omega <- function(alpha, individual = FALSE) {
  alpha <- as.matrix(alpha)
  S <- nrow(alpha)
  # omega function
  omega <- function(S, Sigma) {
    m <- matrix(0, S, 1)
    a <- matrix(0, S, 1)
    b <- matrix(Inf, S, 1)
    d <- pmvnorm(lower = rep(0, S), upper = rep(Inf, S), mean = rep(0, S), sigma = Sigma)
    if(individual){
      out <- d[1]^(1/S) # species level
    }else{
      out <- d[1] # community level
    }
    return(out)
  }
  # rule out errors
  f <- function(m) class(try(solve(t(m) %*% m), silent = T)) == "matrix"
  if (all(f(alpha) == FALSE)) {
    return(0)
  }
  else {
    Sigma <- solve(t(alpha) %*% alpha)
    return(omega(S, Sigma))
  }
}

#' function that normalizes a vector in the L2 norm
#' @param a the original vector
#' @param norm specify L^n norm
#' @return the normalized vector
#' @export
normalize <- function(a, norm = 2) {
  a / sqrt(sum(a^2))
}

#' function that normalizes the spanning vectors of the feasibility domain in the L2 norm
#' @param alpha interaction matrix, diag(alpha)=-1
#' @return the normalized spanning vectors
#' @export
generate_span_vectors <- function(alpha) {
  apply(-alpha, 2, function(x) normalize(x))
}

#### parmeterize competitive interaction matrix with metabolic theory
interaction_matrix_random <- function(M,num,alpha,E){
  Inte <- runif(num*num)
  Inte <- matrix(Inte, ncol = num, nrow = num)
  for (a in 1:length(M)) {
    for (b in 1:length(M)) {
      Inte[a,b] <- (M[b]/M[a])^(alpha)*exp((E[a]-E[b]))*(1/sqrt(num))
    }
  }
  diag(Inte) <- 1 # per niche framework
  Inte <- abs(Inte)*-1 #competition
  return(Inte)
}


#INITIALIZATION
num = 5 # number of species
loop <- 10^5 ## ITERATIONS
### range of alphas to study
ralpha <- c(-2,-1.75,-1.5,-1.25,-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1,1.25,1.5,1.75,2)
prob <- rep(0,length(ralpha))
for (i in 1:length(ralpha)) {
  size <- 0
  alpha <- ralpha[i]
  for (j in 1:loop) {
    M <- runif(num,0,1) ## MASS
    M <- M / sum(M)
    E <- runif(num,0.5,0.7) # activation energy
    A <- interaction_matrix_random(M,num,alpha,E)
    ss <- (2^num*Omega(A))^(1/num) ### prob. of feasibility per population
    size <- size + ss # counter
  }
  size <- size / loop # average
  prob[i] <- size #save results
}
plot(ralpha,prob,xlab ='Alpha',ylab = 'Probability',type='o',col='darkorange',lwd=2,cex.lab=1,axis.ticks=3)
axis(1, at = seq(-2,2,by = 0.25),labels = c(-2,-1.75,-1.5,-1.25,-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1,1.25,1.5,1.75,2))
