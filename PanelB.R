rm(list=ls()) ## clean workspace
set.seed(1) ## for reproduction

######## FUNCTIONS ##################


#### parmeterize competitive interaction matrix with metabolic theory
interaction_matrix_random <- function(M,num,alpha,E){
  Inte <- runif(num*num)
  Inte <- matrix(Inte, ncol = num, nrow = num)
  for (a in 1:length(M)) {
    for (b in 1:length(M)) {
      Inte[a,b] <- (M[b]/M[a])^(alpha)*exp((E[a]-E[b]))*(1/sqrt(num))
    }
  }
  diag(Inte) <- 1  # per niche framework
  Inte <- abs(Inte)*-1  #competition
  return(Inte)
}


### FUNCTION TO CALCULATE THE DISTANCE
theta <- function(H,r){
  r_c <- v_centroid(H)
  out <- acos(sum(r_c*r)/(sqrt(sum(r^2))*sqrt(sum(r_c^2))))*180/pi
  return(out)
}


### FUNCTION TO GET THE CENTROID
v_centroid <- function(G){
  n <- nrow(G)
  D <- diag(1/sqrt(diag(t(G)%*%G)))
  G_n <- G %*% D
  v_c <- rowSums(G_n) /n 
  v_c <- t(t(v_c))
  return(v_c)
}

### Function to parameterize vector K
v_mt <- function(beta,num,M){
  v <- M^(beta)
  v <- v / sqrt(sum(v^2)) #normalization on the sphere
  v
}

#INITIALIZATION
num = 5 # number of species
loop <- 10^5 ## ITERATIONS
### range of alphas to study
ralpha <- c(-2,-1.75,-1.5,-1.25,-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1,1.25,1.5,1.75,2)
angle <- rep(0,length(ralpha))

## generate matrices abd analysis
for (i in 1:length(ralpha)) {
  alpha <- ralpha[i]
  for (j in 1:loop) {
    M <- runif(num,0,1) ## MASS
    M <- M / sum(M)
    E <- runif(num,0.5,0.7) ## ACTIVATION ENERGY
    A <- interaction_matrix_random(M,num,alpha,E) ## changing ALPHA
    r_mt <- v_mt(beta=-3/4,num,M) ## specify the beta to study
  }
  angle[i] <- theta(-A,r_mt)
}


#### PLOT
plot(ralpha,angle,xlab ='Alpha',xaxt = "n",ylab = 'Average distance',type='o',col='dodgerblue',cex.lab=1,lwd=1.5,axis.ticks=3,ylim = c(0,90))
axis(1, at = seq(-2,2,by = 0.25),labels = c(-2,-1.75,-1.5,-1.25,-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1,1.25,1.5,1.75,2))
