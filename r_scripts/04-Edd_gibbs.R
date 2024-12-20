
library(coda)

############ 2.3 Gibbs algorithm 
############ 
gibbs_sampler <- function(N = 50000, burnin = 20000, my_seed = 2025,
                          theta1_init =0.1,theta2_init = 0.1){
  #set seed
  set.seed(my_seed)
  #sample storage
  theta1_sample <- theta2_sample  <- numeric(N)
  #initialise values
  theta1_sample[1] = theta1_init
  theta2_sample[1] = theta2_init
  
  for (i in 2:N) {
    theta1_sample[i] <- rnorm(1,
                     mean = -(cos(2 * pi + 0.3) * theta2_sample[i - 1] + 0.3) /
                       (2 * (8 * theta2_sample[i - 1]^2 + 0.5)),
                     sd = 1 / sqrt(2 * (8 * theta2_sample[i - 1]^2 + 0.5)))
  
    theta2_sample[i] <- rnorm(1,
                     mean = -(cos(2 * pi + 0.3) * theta1_sample[i] + 0.2) /
                       (2 * (8 * theta1_sample[i]^2 + 0.5)),
                     sd = 1 / sqrt(2 * (8 * theta1_sample[i]^2 + 0.5)))
}

# Return samples after burn-in
  theta1_sample = theta1_sample[(burnin + 1):N]
  theta2_sample = theta2_sample[(burnin + 1):N]
# set an out list
  outlist <- list(theta1 = theta1_sample,
                  theta2 = theta2_sample,
                  N = N,
                  burnin = burnin)
  return(outlist)
}


# Summary and assessing convergence 
MCMCsample <- gibbs_sampler()
theta1 <- MCMCsample$theta1
theta2 <- MCMCsample$theta2

#Histograms
par(mfrow=c(1,2))
hist(theta1, breaks = seq(min(theta1), max(theta1), length.out = 200), freq = FALSE,
     xlim = c(-5,5), ylim=c(0,1), ylab=expression(p(theta[1])), xlab=expression(theta[1]), 
     col="green", main="Gibbs algorithm (N=30000 samples)", cex.main=0.8)


hist(theta2, breaks = seq(min(theta2), max(theta2), length.out = 200), freq = FALSE,
     xlim = c(-5,5), ylim=c(0,1), ylab=expression(p(theta[2])), xlab=expression(theta[2]),
     col="green", main="Gibbs algorithm (N=30000 samples)", cex.main=0.8)


####Trace Plots
par(mfrow = c(1,2))
coda::traceplot(as.mcmc(theta1), main = "theta1", col = "orange")
coda::traceplot(as.mcmc(theta2), main = "theta2", col = "orange")

###summary of 
summary(theta1)
HPDinterval(as.mcmc(theta1))
quantile(theta1, probs = c(0.025, 0.975))

summary(theta2)
HPDinterval(as.mcmc(theta2))
quantile(theta1, probs = c(0.025, 0.975))

## Question 5 ---------------------
#--- Geweke diagnostic
coda::geweke.diag(as.mcmc(theta1))
coda::geweke.diag(as.mcmc(theta2))

## Question 6 -------------
mean(theta1)
quantile(theta1, probs = c(0.025, 0.975))

## Question 7 -------------
#estimate the probability P(theta1 > 0.5)
(1+sum(theta1>0.5))/(30000+1)


