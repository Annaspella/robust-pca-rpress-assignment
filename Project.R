#Model
library(rospca)

# Generate training data with contaminated matrix

#gen.data(coeff, n = 100, eps = 0.1, sig = 3, snr = 1/20, seed = 837)
#picontam <- 0.10 # Percentage of contaminated observations


model <- function(n, mu =0, sigma, pi = 0.10, mu_tilde)
{
  x <- matrix(NA, n, p)
  v <- rbinom(n, 1, pi)
  for (i in 1:n){
    if (v[i]==0)
      x[i,] <- rnorm(p, mu, sqrt(sigma[i,i]))
    else
      x[i,] <- rnorm(p, mu_tilde, sqrt(sigma[i,i]))
  }
  
  return(x) 
}

# Try to understand how to select mu_tilde and sigma

### Models with contamination
#Setting (1)
set.seed(0)
n <- 30 # Number of observations
p <- 50  # Number of dimensions

x1_out <- array(NA, dim=c(n,p,10))

for(j in 1:10)
  x1_out[,,j] <- model(n=n, sigma=diag(rep(3,p)), mu_tilde=c(rep(0,5), 50, rep(0,(p-6))))

dim(x1_out)
summary(x1_out)

#setting (2)
set.seed(666)
n <- 100
p <- 500

x2_out <- array(NA, dim=c(n,p,10))

for(j in 1:10) 
  x2_out[,,j] <- model(n=n, sigma=diag(rep(3,p)), mu_tilde=c(rep(0,5), 50, rep(0,(p-6))))

dim(x2_out)
summary(x2_out)

#Models without contamination (pi=0)
#Setting (1)
set.seed(111)
n <- 30 # Number of observations
p <- 50  # Number of dimensions

x1 <- array(NA, dim=c(n,p,10))

for(j in 1:10)
  x1[,,j] <- model(n=n, sigma=diag(rep(3,p)), pi=0, mu_tilde=c(rep(0,5), 50, rep(0,(p-6))))

dim(x1)
summary(x1)

#setting (2)
set.seed(111)
n <- 100
p <- 500

x2 <- array(NA, dim=c(n,p,10))

for(j in 1:10)
  x2[,,j] <- model(n=n, sigma=diag(rep(3,p)), pi=0, mu_tilde=c(rep(0,5), 50, rep(0,(p-6))))

dim(x2)
summary(x2)


#### Model with datagen:

ppp <- dataGen(m = 10, n = 100, p = 500, a = c(0.9,0.5,0), bLength = 4, SD = c(10,5,2),
        eps = 0, seed = TRUE)


#### Calculating R-PRESS


r1 <- robpca (x1_out[,,1], k = 0, kmax = 10, alpha = 0.75, h = NULL, mcd = FALSE,
      ndir = "all", skew = FALSE)



