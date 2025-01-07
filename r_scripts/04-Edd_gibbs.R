
library(coda)

#####set working directory
library(here)
setwd(paste0(here(),"/bayesian-data-analysis"))
############ 2.3 Gibbs algorithm 
############ 
gibbs_sampler <- function(N = 50000, burnin = 20000, my_seed = 2025,
                          theta1_init = 0.1, theta2_init = 0.1){
  #set seed
  set.seed(my_seed)
  #sample storage
  theta1_sample <- theta2_sample  <- numeric(N)
  #initialise values
  theta1_sample[1] = theta1_init
  theta2_sample[1] = theta2_init
  
  for (i in 2:N) {
    theta1_sample[i] <- rnorm(1,
                     mean = (cos(2 * pi + 0.3) * theta2_sample[i - 1] + 0.3) /
                       ((16 * theta2_sample[i - 1]^2) + 1),
                     sd = 1 / sqrt((16 * theta2_sample[i - 1]^2) + 1))
  
    theta2_sample[i] <- rnorm(1,
                     mean = (cos(2 * pi + 0.3) * theta1_sample[i] + 0.2) /
                       ((16*theta1_sample[i]^2)+ 1),
                     sd = 1 / sqrt((16 * theta1_sample[i]^2) + 1))
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



####Trace Plots
png("pictures/fig01-gs-traceplot.png", width = 18, 
    height = 9, units = "cm", res = 300)
par(mfrow = c(1,2))
coda::traceplot(as.mcmc(theta1), main = "Theta 1", col = "#51127c")
coda::traceplot(as.mcmc(theta2), main = "Theta 2", col = "#fc8961")
dev.off()

#######Histograms
png("pictures/fig02-gs-histogram.png", width = 18, 
    height = 9, units = "cm", res = 300)
par(mfrow=c(1,2))
hist(theta1, breaks = seq(min(theta1), max(theta1), length.out = 200), freq = FALSE,
     xlim = c(-5,5), ylim=c(0,1), ylab=expression(p(theta[1])), xlab=expression(theta[1]), 
     col="#51127c", cex.main=0.8)

hist(theta2, breaks = seq(min(theta2), max(theta2), length.out = 200), freq = FALSE,
     xlim = c(-5,5), ylim=c(0,1), ylab=expression(p(theta[2])), xlab=expression(theta[2]),
     col="#fc8961", cex.main=0.8)
dev.off()


# ----- Q-Q plots (first half and second half of theta 1 and theta 2)
png("pictures/fig03-gs-qqplots.png", width = 18, 
    height = 9, units = "cm", res = 300)

par(mfrow = c(1, 2))

theta1_1st <- theta1[1:(50000 - 20000)/2]
theta1_2nd <- theta1[((50000 - 20000)/2 + 1):(50000 - 20000)]

qqplot(theta1_1st, theta1_2nd, xlab = "First half of Theta 1",
       ylab = "Second half of Theta 1", col = "#51127c")

theta2_1st <- theta2[1:(50000 - 20000)/2]
theta2_2nd <- theta2[((50000 - 20000)/2 + 1):(50000 - 20000)]
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

runmean_theta1 <- runmean(theta1)
runmean_theta2 <- runmean(theta2)

png("pictures/fig04-gs-runmean.png", width = 18, 
    height = 9, units = "cm", res = 300)

par(mfrow = c(1, 2))

plot(seq_len(50000 - 20000), runmean_theta1, 
     type = "l", xlab = "Iterations", 
     ylab = "Theta 1", col = "#51127c", lwd = 2, yaxt = "n")
axis(2, at = seq(-0.8, 0.2, by = 0.05), las = 2)

plot(seq_len(50000 - 20000), runmean_theta2, 
     type = "l", xlab = "Iterations", 
     ylab = "Theta 2", col = "#fc8961", lwd = 2, yaxt = "n")
axis(2, at = seq(-0.05, 0.25, by = 0.05), las = 2)

dev.off()


 #of 
summary(theta1)
HPDinterval(as.mcmc(theta1))
quantile(theta1, probs = c(0.025, 0.975))

summary(theta2)
HPDinterval(as.mcmc(theta2))
quantile(theta1, probs = c(0.025, 0.975))

## Question 5 ---------------------
#--------Geweke diagnostic
 coda::geweke.diag(as.mcmc(theta1))
 coda::geweke.diag(as.mcmc(theta2))

###-------geweke diag plot
t1 <- as.data.frame(theta1)
 t2 <- as.data.frame(theta2)
 thetadf <- cbind(t1,t2)
 
png("pictures/fig05-gs-geweke.png", width = 18, 
height = 9, units = "cm", res = 300)
coda::geweke.plot(as.mcmc((thetadf)))

dev.off()


## Question 6 -------------
mean(theta1)
quantile(theta1, probs = c(0.025, 0.975))

## Question 7 -------------
#estimate the probability P(theta1 > 0.5)
(1+sum(theta1>0.5))/(30000+1)


