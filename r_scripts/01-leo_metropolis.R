## load libraries -------------------------------------------------------------
library(mvtnorm)
library(coda)
library(ggplot2)
library(HDInterval)

set.seed(1993)

## Question 2 - Newton-Raphson algorithm -------------------------------------

# define constants
A <- 16
B <- cos(2 * pi + 0.3)
C1 <- 0.3
C2 <- 0.2

# Define the gradient function in the global environment
gradient_theta <- function(theta = c(1, 2)) {
  t1 <- theta[1]
  t2 <- theta[2]
  grad1 <- -A * t1 * t2^2 - t1 + B * t2 + C1
  grad2 <- -A * t1^2 * t2 - t2 + B * t1 + C2
  return(c(grad1, grad2))
}

# Define the Hessian function in the global environment
hessian_theta <- function(theta = c(1, 2)) {
  t1 <- theta[1]
  t2 <- theta[2]
  h11 <- -A * t2^2 - 1
  h12 <- h21 <- -2 * A * t1 * t2 + B
  h22 <- -A * t1^2 - 1
  return(matrix(c(h11, h12, h21, h22), nrow = 2))
}

NewtonRaphson <- function(theta0 = c(1, 2), tolerance = 1e-10, iteration = 1000) {
  for (i in 1:iteration) {
    grad_theta0 <- gradient_theta(theta0)
    hes_theta0 <- hessian_theta(theta0)
    theta1 <- theta0 - solve(hes_theta0, grad_theta0) # Calculate next value of theta
    distance <- sum(abs(theta1 - theta0))
    if (distance < tolerance) { # Once distance between theta1 and theta0 gets sufficiently small, output result
      result <- round(theta1, 5)
      return(list("mode theta" = result, paste("Converged at iteration", i)))
    }
    # If convergence not reached, set theta1 as theta0 and continue
    theta0 <- theta1
  }
  print('Newton-Raphson algorithm did not converge. Choose a larger number of iterations.')
}

theta_NR_star <- NewtonRaphson()$'mode theta'
theta_NR_star


## Question 3: Calculate covariance matrix of Laplace. ----------------------
# It is the inverse of the negative Hessian matrix
sigma_star <- round(solve(-hessian_theta(theta = theta_NR_star)),5)
sigma_star

# Calculate the mean
grad_theta_nr <- gradient_theta(theta_NR_star)
hes_theta_nr <- hessian_theta(theta_NR_star)
mean <- round(theta_NR_star - solve(hes_theta_nr,grad_theta_nr),5)
mean


## Question 4: Random walk Metropolis algorithm.------------------------------ 

Metropolis <- function(M = 51000, burn_in = 1000, seed = 1993, 
                       theta_start = c(0.29970, 0.19955), c_tune = 2 )
{
  set.seed(seed)
  
  # Joint log posterior distribution
  log_posterior <- function(theta)
  {
    t1 <- theta[1]
    t2 <- theta[2]
    return(-1/2*(A*t1^2*t2^2 + t1^2 + t2^2 - 2*B*t1*t2 - 2*C1*t1 - 2*C2*t2))
  }
  
  # Starting value for theta 1 and theta 2
  theta <- array(dim = c(2, M))
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
    theta_prop <- theta[,i - 1] + rmvnorm(n = 1, 
                                        mean = c(0,0), sigma = sigma_tilde)
    # Log ratio for accept-reject decision
    logr <- log_posterior(theta_prop) - log_posterior(theta[,i - 1])
    
    # Accept or reject
    u <- runif(1)
    if (logr >= 0 || u <= exp(logr))
    { #Accept
      theta[,i] <- theta_prop
      n_accept <- n_accept + 1
    } 
    else 
    { #Reject, stay where theta is
      theta[,i] <- theta[, i - 1]
    }
  }
  # Exclude burn-in
  theta <- theta[,-c(1:burn_in)]
  # Output
  accept_rate <- round(n_accept/(M - 1), digits = 3) * 100
  rownames(theta) <- c("Theta1","Theta2")
  output <- list(theta = theta,
                 M = M, 
                 burn_in = burn_in,
                 log_posterior = log_posterior)
  print(paste("Acceptance rate:", accept_rate,"%"))
  return(invisible(output))
}

MCMC_run <- Metropolis(c_tune = 2.8)
burn_in <- MCMC_run$burn_in
theta <- MCMC_run$theta
M <- MCMC_run$M

### summary stats -----------------------------------------------------------

# Summary stats of theta 1 and theta 2
summary(as.mcmc(theta[1,]))
# Mean             SD       Naive SE Time-series SE 
# 0.165480       0.685817       0.003067       0.011127 
# 2.5%      25%      50%      75%    97.5% 
# -1.10059 -0.21200  0.08443  0.50623  1.72628
hdi(as.mcmc(theta[1,]))
# var1
# lower -1.100592
# upper  1.723399
# attr(,"credMass")
# [1] 0.95

summary(as.mcmc(theta[2,]))
# Mean             SD       Naive SE Time-series SE 
# 0.100503       0.647716       0.002897       0.010574 
# 2.5%      25%      50%      75%    97.5% 
# -1.21044 -0.24057  0.06609  0.41506  1.57926
hdi(as.mcmc(theta[2,]))
# var1
# lower -1.231357
# upper  1.547863
# attr(,"credMass")
# [1] 0.95

### -------Assess convergence-------------------------------------------------

# ------Trace plots of theta 1 and theta 2
png("pictures/fig01-mh-traceplot.png", width = 18, 
    height = 9, units = "cm", res = 300)

par(mfrow = c(1, 2))
coda::traceplot(as.mcmc(theta[1,]), main = "Theta 1", col = "#51127c")
coda::traceplot(as.mcmc(theta[2,]), main = "Theta 2", col = "#fc8961")

dev.off()

# -----Autocorrelation plots
png("pictures/fig02-mh-autocorr.png", width = 18, 
    height = 9, units = "cm", res = 300)

par(mfrow = c(1, 2))
coda::autocorr.plot(as.mcmc(theta[1,]), lag.max = 40, 
                    main = "Theta 1", lwd = 2, col = "#51127c",
                    auto.layout = F)
coda::autocorr.plot(as.mcmc(theta[2,]), lag.max = 40, 
                    main = "Theta 2", lwd = 2, col = "#fc8961",
                    auto.layout = F)

dev.off()

# ----- Q-Q plots (first half and second half of theta 1 and theta 2)
png("pictures/fig03-mh-qqplots.png", width = 18, 
    height = 9, units = "cm", res = 300)

par(mfrow = c(1, 2))

theta1_1st <- theta[1, 1:(M - burn_in)/2]
theta1_2nd <- theta[1, ((M - burn_in)/2 + 1):(M - burn_in)]
qqplot(theta1_1st, theta1_2nd, xlab = "First half of Theta 1",
       ylab = "Second half of Theta 1", col = "#51127c")

theta2_1st <- theta[2, 1:(M - burn_in)/2]
theta2_2nd <- theta[2, ((M - burn_in)/2 + 1):(M - burn_in)]
qqplot(theta2_1st, theta2_2nd, xlab = "First half of Theta 2",
       ylab = "Second half of Theta 2", col = "#fc8961")

dev.off()

# ------ Running mean to assess stability of the mean
runmean <- function(x){
  nlen <- length(x)
  mean_iteration <- c()
  mean_iteration[1] <- x[1]
  for (j in 2:nlen) {
    mean_iteration[j] <- mean(x[1:j])
  }
  return(mean_iteration)
}

runmean_theta1 <- runmean(theta[1,])
runmean_theta2 <- runmean(theta[2,])

png("pictures/fig04-mh-runmean.png", width = 18, 
    height = 9, units = "cm", res = 300)

par(mfrow = c(1, 2))

plot(seq_len(M - burn_in), runmean_theta1, 
     type = "l", xlab = "Iterations", 
     ylab = "Theta 1", col = "#51127c", lwd = 2, yaxt = "n")
axis(2, at = seq(-0.8, 0.2, by = 0.05), las = 2)

plot(seq_len(M - burn_in), runmean_theta2, 
     type = "l", xlab = "Iterations", 
     ylab = "Theta 2", col = "#fc8961", lwd = 2, yaxt = "n")
axis(2, at = seq(-0.05, 0.25, by = 0.05), las = 2)

dev.off()

# ----- Geweke diagnostic
png("pictures/fig05-mh-geweke.png", width = 18, 
    height = 9, units = "cm", res = 300)

coda::geweke.diag(as.mcmc(theta[1,]))
coda::geweke.diag(as.mcmc(theta[2,]))
coda::geweke.plot(as.mcmc(t(theta)))

dev.off()

# ----- Gelman_Rubin_Diagnostic
# Run some other chains
MCMC_run1 <- Metropolis(c_tune = 2.8, 
                        theta_start = c(0.3, 0.2))
MCMC_run2 <- Metropolis(c_tune = 2.8, 
                        theta_start = c(0.2, 0.1))
MCMC_run3 <- Metropolis(c_tune = 2.8, 
                        theta_start = c(0.4, 0.3))

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

# export gelman plots
png("pictures/fig06-mh-gelman.png", width = 18, 
    height = 9, units = "cm", res = 300)

par(mfrow = c(1, 2))
coda::gelman.plot(theta1_obj, main = "Theta 1", auto.layout = F)
coda::gelman.plot(theta2_obj, main = "Theta 2", auto.layout = F)

dev.off()

## ---- Question 5: P(theta1/theta2) >0.45 ------------------

ratio_theta12 <- theta_run1[1,]/theta_run1[2,]
p_ratio_045 <- sum(ifelse(ratio_theta12 > 0.45, 1, 0))/length(ratio_theta12)
p_ratio_045