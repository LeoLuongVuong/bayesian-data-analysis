# Required library
library(mvtnorm)
library(coda)

###  ------------------------------------- PART 1 -------------------------------------
### Part 1.1 - Gibbs algorithm -------------------------------------

## Question 4 ----------------
Gibbs <- function(M = 35000, burn = 5000, seed = 2025){
  # Data
  y <- c(4, 1, 5, 14, 3, 19, 7, 6)
  t <- c(95, 16, 63, 126, 6, 32, 16, 19)
  # Constant
  alpha <- 1.8
  gamma <- 0.01
  delta <- 1
  # Gibbs sampling
  set.seed(seed)
  beta_sample <-lambda1_sample <- lambda2_sample <- lambda3_sample <- lambda4_sample <- lambda5_sample <- lambda6_sample <- lambda7_sample <- lambda8_sample <- c()
  lambda_init <- c(rep(0.1,8))
  M <- M
  burnin <- burn
  for(m in 1:M)
  {
    beta_sample[m] <- rgamma(n = 1, shape = 8*alpha + gamma, rate = sum(lambda_init) + delta)
    lambda1_sample[m] <- rgamma(n = 1, shape = y[1] + alpha, rate = t[1] + beta_sample[m])
    lambda2_sample[m] <- rgamma(n = 1, shape = y[2] + alpha, rate = t[2] + beta_sample[m])
    lambda3_sample[m] <- rgamma(n = 1, shape = y[3] + alpha, rate = t[3] + beta_sample[m])
    lambda4_sample[m] <- rgamma(n = 1, shape = y[4] + alpha, rate = t[4] + beta_sample[m])
    lambda5_sample[m] <- rgamma(n = 1, shape = y[5] + alpha, rate = t[5] + beta_sample[m])
    lambda6_sample[m] <- rgamma(n = 1, shape = y[6] + alpha, rate = t[6] + beta_sample[m])
    lambda7_sample[m] <- rgamma(n = 1, shape = y[7] + alpha, rate = t[7] + beta_sample[m])
    lambda8_sample[m] <- rgamma(n = 1, shape = y[8] + alpha, rate = t[8] + beta_sample[m])
    lambda_init[1] <- lambda1_sample[m]
    lambda_init[2] <- lambda2_sample[m]
    lambda_init[3] <- lambda3_sample[m]
    lambda_init[4] <- lambda4_sample[m]
    lambda_init[5] <- lambda5_sample[m]
    lambda_init[6] <- lambda6_sample[m]
    lambda_init[7] <- lambda7_sample[m]
    lambda_init[8] <- lambda8_sample[m]
  }
  # Discard burn in
  beta_sample <- beta_sample[-c(1:burnin)]
  lambda1_sample <- lambda1_sample[-c(1:burnin)]
  lambda2_sample <- lambda2_sample[-c(1:burnin)]
  lambda3_sample <- lambda3_sample[-c(1:burnin)]
  lambda4_sample <- lambda4_sample[-c(1:burnin)]
  lambda5_sample <- lambda5_sample[-c(1:burnin)]
  lambda6_sample <- lambda6_sample[-c(1:burnin)]
  lambda7_sample <- lambda7_sample[-c(1:burnin)]
  lambda8_sample <- lambda8_sample[-c(1:burnin)]
  
  outlist <- list(beta=beta_sample,
                  lambda1=lambda1_sample,
                  lambda2=lambda2_sample,
                  lambda3=lambda3_sample,
                  lambda4=lambda4_sample,
                  lambda5=lambda5_sample,
                  lambda6=lambda6_sample,
                  lambda7=lambda7_sample,
                  lambda8=lambda8_sample,
                  M=M,
                  burnin = burnin)
  return(outlist)
}

# Summary and assessing convergence 
MCMCrun <- Gibbs()
beta <- MCMCrun$beta
lambda1 <- MCMCrun$lambda1
lambda2 <- MCMCrun$lambda2
lambda3 <- MCMCrun$lambda3
lambda4 <- MCMCrun$lambda4
lambda5 <- MCMCrun$lambda5
lambda6 <- MCMCrun$lambda6
lambda7 <- MCMCrun$lambda7
lambda8 <- MCMCrun$lambda8
M <- MCMCrun$M
burnin <- MCMCrun$burnin

par(mfrow=c(3,3))
hist(beta_sample, breaks = seq(min(beta), max(beta), length.out = 200), freq = FALSE, xlim = c(0,10), ylim=c(0,0.5), ylab=expression(p(beta)), xlab=expression(beta), col="green", main="Gibbs algorithm (M=35000 samples)", cex.main=0.8)
hist(lambda1, breaks = seq(min(lambda1), max(lambda1), length.out = 200), freq = FALSE, xlim = c(0,0.25), ylim=c(0,20), ylab=expression(p(lambda1)), xlab=expression(lambda1), col="green", main="Gibbs algorithm (M=35000 samples)", cex.main=0.8)
hist(lambda2, breaks = seq(min(lambda2), max(lambda2), length.out = 200), freq = FALSE, xlim = c(0,0.5), ylim=c(0,8), ylab=expression(p(lambda2)), xlab=expression(lambda2), col="green", main="Gibbs algorithm (M=35000 samples)", cex.main=0.8)
hist(lambda3, breaks = seq(min(lambda3), max(lambda3), length.out = 200), freq = FALSE, xlim = c(0,0.25), ylim=c(0,15), ylab=expression(p(lambda3)), xlab=expression(lambda3), col="blue", main="Gibbs algorithm (M=35000 samples)", cex.main=0.8)
hist(lambda4, breaks = seq(min(lambda4), max(lambda4), length.out = 200), freq = FALSE, xlim = c(0,0.3), ylim=c(0,15), ylab=expression(p(lambda4)), xlab=expression(lambda4), col="blue", main="Gibbs algorithm (M=35000 samples)", cex.main=0.8)
hist(lambda5, breaks = seq(min(lambda5), max(lambda5), length.out = 200), freq = FALSE, xlim = c(0,1.5), ylim=c(0,3), ylab=expression(p(lambda5)), xlab=expression(lambda5), col="blue", main="Gibbs algorithm (M=35000 samples)", cex.main=0.8)
hist(lambda6, breaks = seq(min(lambda6), max(lambda6), length.out = 200), freq = FALSE, xlim = c(0.2,1.0), ylim=c(0,5), ylab=expression(p(lambda6)), xlab=expression(lambda6), col="red", main="Gibbs algorithm (M=35000 samples)", cex.main=0.8)
hist(lambda7, breaks = seq(min(lambda7), max(lambda7), length.out = 200), freq = FALSE, xlim = c(0.1,0.8), ylim=c(0,5), ylab=expression(p(lambda7)), xlab=expression(lambda7), col="red", main="Gibbs algorithm (M=35000 samples)", cex.main=0.8)
hist(lambda8, breaks = seq(min(lambda8), max(lambda8), length.out = 200), freq = FALSE, xlim = c(0,1), ylim=c(0,5), ylab=expression(p(lambda8)), xlab=expression(lambda8), col="red", main="Gibbs algorithm (M=35000 samples)", cex.main=0.8)

# Trace plots
par(mfrow = c(3,3))
coda::traceplot(as.mcmc(beta), main = "beta", col = "orange")
coda::traceplot(as.mcmc(lambda1), main = "lambda_1", col = "orange")
coda::traceplot(as.mcmc(lambda2), main = "lambda_2", col = "orange")
coda::traceplot(as.mcmc(lambda3), main = "lambda_3", col = "darkgreen")
coda::traceplot(as.mcmc(lambda4), main = "lambda_4", col = "darkgreen")
coda::traceplot(as.mcmc(lambda5), main = "lambda_5", col = "darkgreen")
coda::traceplot(as.mcmc(lambda6), main = "lambda_6", col = "darkblue")
coda::traceplot(as.mcmc(lambda7), main = "lambda_7", col = "darkblue")
coda::traceplot(as.mcmc(lambda8), main = "lambda_8", col = "darkblue")

summary(beta)
HPDinterval(as.mcmc(beta))
quantile(beta, probs = c(0.025, 0.975))

summary(lambda1)
HPDinterval(as.mcmc(lambda1))
quantile(lambda1, probs = c(0.025, 0.975))

summary(lambda2)
HPDinterval(as.mcmc(lambda2))
quantile(lambda2, probs = c(0.025, 0.975))

summary(lambda3)
HPDinterval(as.mcmc(lambda3))
quantile(lambda3, probs = c(0.025, 0.975))

summary(lambda4)
HPDinterval(as.mcmc(lambda4))
quantile(lambda4, probs = c(0.025, 0.975))

summary(lambda5)
HPDinterval(as.mcmc(lambda5))
quantile(lambda5, probs = c(0.025, 0.975))

summary(lambda6)
HPDinterval(as.mcmc(lambda6))
quantile(lambda6, probs = c(0.025, 0.975))

summary(lambda7)
HPDinterval(as.mcmc(lambda7))
quantile(lambda7, probs = c(0.025, 0.975))

summary(lambda8)
HPDinterval(as.mcmc(lambda8))
quantile(lambda8, probs = c(0.025, 0.975))

## Question 5 ---------------------
#--- Geweke diagnostic
coda::geweke.diag(as.mcmc(beta))
coda::geweke.diag(as.mcmc(lambda1))
coda::geweke.diag(as.mcmc(lambda2))
coda::geweke.diag(as.mcmc(lambda3))
coda::geweke.diag(as.mcmc(lambda4))
coda::geweke.diag(as.mcmc(lambda5))
coda::geweke.diag(as.mcmc(lambda6))
coda::geweke.diag(as.mcmc(lambda7))
coda::geweke.diag(as.mcmc(lambda8))

par(mfrow=c(1,3))
theta <- array(dim = c(9, M-burnin))
rownames(theta) <- c("beta","lambda1","lambda2","lambda3","lambda4","lambda5","lambda6","lambda7","lambda8")
theta[1,] <- beta
theta[2,] <- lambda1
theta[3,] <- lambda2
theta[4,] <- lambda3
theta[5,] <- lambda4
theta[6,] <- lambda5
theta[7,] <- lambda6
theta[8,] <- lambda7
theta[9,] <- lambda8
coda::geweke.plot(as.mcmc(t(theta)))

# Heidelberger-Welch stationary diagnostic
heidel.diag(beta, eps=0.1, pvalue=0.05)
heidel.diag(lambda1, eps=0.1, pvalue=0.05)
heidel.diag(lambda2, eps=0.1, pvalue=0.05)
heidel.diag(lambda3, eps=0.1, pvalue=0.05)
heidel.diag(lambda4, eps=0.1, pvalue=0.05)
heidel.diag(lambda5, eps=0.1, pvalue=0.05)
heidel.diag(lambda6, eps=0.1, pvalue=0.05)
heidel.diag(lambda7, eps=0.1, pvalue=0.05)
heidel.diag(lambda8, eps=0.1, pvalue=0.05)

## Question 6 -------------
t[6]*mean(lambda6)
t[6]*quantile(lambda6, probs = c(0.025, 0.975))

## Question 7 -------------
#estimate the probability P(lambda6 > 0.53)
#p monte carlo
(1+sum(lambda6>0.53))/(M-burnin+1)

### Part 1.2 - Metropolis algorithm ----------------------------------------------

## Question 2 - Newton-Raphson algorithm -------------------------------------

# Define constants
A <- 1.9
B <- 0.54
C1 <- 0.5
C2 <- -1.4

#Calculate gradients
gradient_theta <- function(theta=c(1,2))
{
  t1 <- theta[1]
  t2 <- theta[2]
  A <- 1.9
  B <- 0.54
  C1 <- 0.5
  C2 <- -1.4
  grad1 <- -A*t1*t2^2 -t1 + B*t2 + C1
  grad2 <- -A*t1^2*t2 -t2 + B*t1 + C2
  return(c(grad1, grad2))
}

# Check if the gradient calculation is correct
theta <- c(1,2)

f <- expression(-1/2*(A*t1^2*t2^2 + t1^2 + t2^2 - 2*B*t1*t2 - 2*C1*t1 -2*C2*t2))
eval(D(f, 't1'), envir=list(t1=1, t2=2))
eval(D(f, 't2'), envir=list(t1=1, t2=2))
gradient_theta(theta)

#Calculate Hessian matrix
hessian_theta <- function(theta=c(1,2))
{
  t1 <- theta[1]
  t2 <- theta[2]
  A <- 1.9
  B <- 0.54
  h11 <- -A*t2^2 - 1
  h12 <- h21 <- -2*A*t1*t2 + B
  h22 <- -A*t1^2 - 1
  return(matrix(c(h11,h12, h21, h22), nrow=2))
}

# Check Hessian calculation
eval(D(D(f, 't1'),'t1'), envir=list(t1=1, t2=2))
eval(D(D(f, 't1'),'t2'), envir=list(t1=1, t2=2))
hessian_theta()

# Function for Newton-Raphson algorithm
NewtonRaphson2 <- function(theta0=c(1,2), tolerance = 1e-10, iteration = 1000){
  for (i in 1: iteration) {
    grad_theta0 <- gradient_theta(theta0)
    hes_theta0 <- hessian_theta(theta0)
    theta1 <- theta0 - solve(hes_theta0,grad_theta0) #calculate next value of theta
    distance <- sum(abs(theta1 - theta0))
    if (distance < tolerance) { #Once distance between theta1 and theta0 gets sufficiently small, output result
      result <- round(theta1, 5)
      return(list("mode theta" = result, paste("converged at iteration", i)))
    }
    # If convergence not reached set theta1 as theta0 and continue
    theta0 <- theta1
  }
print('Newton-Raphson algorithm did not converge. Choose larger number of iteration')
}

# Function for Newton-Raphson algorithm, version stand-alone
NewtonRaphson <- function(theta0=c(1,2), tolerance = 1e-10, iteration = 1000){
  # Define constants
  A <- 1.9
  B <- 0.54
  C1 <- 0.5
  C2 <- -1.4  
  #Calculate gradients
  gradient_theta <- function(theta=c(1,2))
  {
    t1 <- theta[1]
    t2 <- theta[2]
    grad1 <- -A*t1*t2^2 -t1 + B*t2 + C1
    grad2 <- -A*t1^2*t2 -t2 + B*t1 + C2
    return(c(grad1, grad2))
  }
  #Calculate Hessian matrix
  hessian_theta <- function(theta=c(1,2))
  {
    t1 <- theta[1]
    t2 <- theta[2]
    h11 <- -A*t2^2 - 1
    h12 <- h21 <- -2*A*t1*t2 + B
    h22 <- -A*t1^2 - 1
    return(matrix(c(h11,h12, h21, h22), nrow=2))
  }
  for (i in 1: iteration) {
    grad_theta0 <- gradient_theta(theta0)
    hes_theta0 <- hessian_theta(theta0)
    theta1 <- theta0 - solve(hes_theta0,grad_theta0) #calculate next value of theta
    distance <- sum(abs(theta1 - theta0))
    if (distance < tolerance) { #Once distance between theta1 and theta0 gets sufficiently small, output result
      result <- round(theta1, 5)
      return(list("mode theta" = result, paste("Converged at iteration", i)))
    }
    # If convergence not reached set theta1 as theta0 and continue
    theta0 <- theta1
  }
  print('Newton-Raphson algorithm did not converge. Choose larger number of iteration')
}

theta_NR_star <- NewtonRaphson()$'mode theta'
theta_NR_star

### Question 3: Calculate covariance matrix of Laplace. ----------------------
# It is the inverse of the negative Hessian matrix
sigma_star <- round(solve(-hessian_theta(theta = theta_NR_star)),5)
sigma_star

# Calculate the mean
grad_theta_nr <- gradient_theta(theta_NR_star)
hes_theta_nr <- hessian_theta(theta_NR_star)
mean <- round(theta_NR_star - solve(hes_theta_nr,grad_theta_nr),5)
mean


### Question 4: Random walk Metropolis algorithm.------------------------------ 

Metropolis <- function(M =51000, burn_in = 1000, seed = 1993, 
                       theta_start = c(-0.0553, -1.4216), c_tune=2 )
{
  set.seed(seed)
  
  # Joint log posterior distribution
  log_posterior <- function(theta)
  {
    t1 <- theta[1]
    t2 <- theta[2]
    return(-1/2*(A*t1^2*t2^2 + t1^2 + t2^2 - 2*B*t1*t2 - 2*C1*t1 -2*C2*t2))
  }
  
  # Starting value for theta 1 and theta 2
  theta <- array(dim=c(2,M))
  theta[, 1] <- theta_start
  n_accept <- 0
  
  # Gaussian proposal
  # Covariance matrix
  theta_NR_star <- NewtonRaphson()$'mode theta' #The NewtonRaphson function built in Q2
  sigma_star <- round(solve(-hessian_theta(theta = theta_NR_star)),5)
  sigma_tilde <- c_tune*sigma_star
  
  # Metropolis loop
  for (i in 2:M)
  {
    # New proposed theta
    theta_prop <- theta[,i-1] + rmvnorm(n = 1, 
                                        mean = c(0,0), sigma = sigma_tilde)
    # Log ratio for accept-reject decision
    logr <- log_posterior(theta_prop) - log_posterior(theta[,i-1])
    
    # Accept or reject
    u <- runif(1)
    if(logr >= 0 || u <= exp(logr))
    { #Accept
      theta[,i] <- theta_prop
      n_accept <- n_accept + 1
    } 
    else 
    { #Reject, stay where theta is
      theta[,i] <- theta[, i-1]
    }
  }
  # Exclude burn-in
  theta <- theta[,-c(1:burn_in)]
  # Output
  accept_rate <- round(n_accept/(M-1), digits = 3) * 100
  rownames(theta) <- c("Theta1","Theta2")
  output <- list(theta= theta,
                 M = M, 
                 burn_in = burn_in,
                 log_posterior = log_posterior)
  print(paste("Acceptance rate:", accept_rate,"%"))
  return(invisible(output))
}

MCMC_run <- Metropolis(c_tune = 6.2, M = 60000, burn_in = 10000)
burn_in <- MCMC_run$burn_in
theta <- MCMC_run$theta
M <- MCMC_run$M

# Summary stats of theta 1 and theta 2
summary(as.mcmc(theta[1,]))
summary(as.mcmc(theta[2,]))


# -------Assess convergence-------------------------------------------------

# ------Trace plots of theta 1 and theta 2
coda::traceplot(as.mcmc(theta[1,]), main = "Theta 1", col = "salmon")
coda::traceplot(as.mcmc(theta[2,]), main = "Theta 2", col = "royalblue")



# -----Autocorrelation plots
coda::autocorr.plot(as.mcmc(theta[1,]), lag.max = 40, 
                    main = "Theta 1", lwd = 3, col = "salmon")
coda::autocorr.plot(as.mcmc(theta[2,]), lag.max = 40, 
                    main = "Theta 2", lwd = 3, col = "royalblue")

# ----- Q-Q plots (first half and second half of theta 1 and theta 2)

theta1_1st <- theta[1, 1:(M-burn_in)/2]
theta1_2nd <- theta[1, ((M-burn_in)/2+1) :(M-burn_in)]
qqplot(theta1_1st, theta1_2nd, xlab = "First half of Theta 1",
       ylab = "Second half of Theta 1", col = "salmon")

theta2_1st <- theta[2, 1:(M-burn_in)/2]
theta2_2nd <- theta[2, ((M-burn_in)/2+1) :(M-burn_in)]
qqplot(theta2_1st, theta2_2nd, xlab = "First half of Theta 2",
       ylab = "Second half of Theta 2", col = "royalblue")

# ------ Running mean to assess stability of the mean
runmean <- function(x){
  nlen <- length(x)
  mean_iteration <- c()
  mean_iteration[1] <- x[1]
  for(j in 2:nlen){
    mean_iteration[j] <- mean(x[1:j])
  }
  return(mean_iteration)
}

runmean_theta1 <- runmean(theta[1,])
runmean_theta2 <- runmean(theta[2,])
plot(seq_len(M-burn_in), runmean_theta1, 
     type = "l", xlab = "Iterations", ylab = "Theta 1", col = "salmon", lwd = 2)
plot(seq_len(M-burn_in), runmean_theta2, 
     type = "l", xlab = "Iterations", ylab = "Theta 2", col = "royalblue", lwd = 2)

# ----- Geweke diagnostic
coda::geweke.diag(as.mcmc(theta[1,]))
coda::geweke.diag(as.mcmc(theta[2,]))
coda::geweke.plot(as.mcmc(t(theta)))

# ----- Gelman_Rubin_Diagnostic
# Run some other chains
MCMC_run1 <- Metropolis(c_tune = 6.2, M = 60000, burn_in = 10000, seed = 1993, 
                         theta_start = c(-0.0553, -1.4216))
MCMC_run2 <- Metropolis(c_tune = 6.2, M = 60000, burn_in = 10000, seed = 1994, 
                        theta_start = c(0.05, 1.4))
MCMC_run3 <- Metropolis(c_tune = 6.2, M = 60000, burn_in = 10000, seed = 1995, 
                        theta_start = c(1, 2))
theta_run1 <- MCMC_run1$theta
theta_run2 <- MCMC_run2$theta
theta_run3 <- MCMC_run3$theta

theta1_obj <- mcmc.list(chain1 = mcmc(theta_run1[1,]), 
                        chain2 = mcmc(theta_run2[1,]),
                        chain3 = mcmc(theta_run1[1,]))

theta2_obj <- mcmc.list(chain1 = mcmc(theta_run1[2,]), 
                        chain2 = mcmc(theta_run2[2,]),
                        chain3 = mcmc(theta_run1[2,]))
coda::gelman.diag(theta1_obj)
coda::gelman.plot(theta1_obj, main = "Theta 1")
coda::gelman.plot(theta2_obj, main = "Theta 2")
 
### ---- Question 5: P(theta1/theta2) >0.45 ------------------

ratio_theta12 <- theta_run1[1,]/ theta_run1[2,]
p_ratio_045 <- sum(ifelse(ratio_theta12 > 0.45, 1, 0))/ length(ratio_theta12)
p_ratio_045


###  ------------------------------------- PART 2 -------------------------------------

###  ---- Part2: Question 3------------------------------------------------------------

# Load libraries and data
library(ggplot2)
library(coda)
library(readr)
library(tidyverse)
library(coda)
setwd("G:/My Drive/Master of Statistics and Data Science/Baysian Data Analysis DL")
GSPS <- readRDS("Project/GSPS.Rdata")

# Load the coda files for the MCMC chains of the predicted values
predict1 <- read.coda("Project/Q3/predict1.txt", "Project/Q3/predict_ind.txt")
predict2 <- read.coda("Project/Q3/predict2.txt", "Project/Q3/predict_ind.txt")
predict3 <- read.coda("Project/Q3/predict3.txt", "Project/Q3/predict_ind.txt")

# Prepare the MCMC object containing multiple chains
mcmcpredict <- mcmc.list(predict1, predict2, predict3)
summary(mcmcpredict)

# Extract the median values
median_values <- apply(as.matrix(mcmcpredict), 2, median)
median_values

# Compare the observed with the predicted values
observed <- GSPS$working
predicted <- median_values
data <- table(predicted,observed)
data

ggplot(df, aes(x = observed, y = predicted, fill = Freq)) +
  geom_tile(color = "white") + geom_text(aes(label = Freq)) +
  scale_fill_gradient(low = "red", high = "darkgreen")
labs(x = "Observed",
     y = "Predicted",
     fill = "Frequency")