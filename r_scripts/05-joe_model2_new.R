library(R2jags)
library(coda)
library(readxl)
library(pacman)
library(MCMCvis)
library(mcmcplots)

DDKdata <- read.csv("bda_data.csv")



datajags2 <- list(
  n_munic = length(unique(DDKdata$naam)),
  n_age = length(unique(DDKdata$age)),
  n_sex = length(unique(DDKdata$sex)),
  Niag = array(c(DDKdata$invited), dim = c(300, 5, 2)),
  Yiag = array(c(DDKdata$participant), dim = c(300, 5, 2)),
  male = array(c(DDKdata$male), dim = c(300, 5, 2)),
  female = array(c(DDKdata$female), dim = c(300, 5, 2)),
  age1 = array(c(DDKdata$age1), dim = c(300, 5, 2)),
  age2 = array(c(DDKdata$age2), dim = c(300, 5, 2)),
  age3 = array(c(DDKdata$age3), dim = c(300, 5, 2)),
  age4 = array(c(DDKdata$age4), dim = c(300, 5, 2)),
  age5 = array(c(DDKdata$age5), dim = c(300, 5, 2))
)

n_munic = length(unique(DDKdata$naam))
n_age = length(unique(DDKdata$age))
n_sex = length(unique(DDKdata$sex))
agec = array(c(DDKdata$agec), dim = c(300, 5, 2))
gender = array(c(DDKdata$sex), dim = c(300, 5, 2))


sink("model2.txt")
cat("
# ...
model {
  for (i in 1:n_munic) {        # Loop over municipalities
    for (j in 1:n_age) {        # Loop over age groups
      for (k in 1:n_sex) {      # Loop over genders

        # binomial likelihood
        # distribution: participation rate
        Yiag[i,j,k] ~ dbin(pi[i,j,k], Niag[i,j,k])  # Likelihood
        # link function:logit of participation rate
        logit(pi[i,j,k]) <- alpha + beta1[j]*age1[i,j,k] + beta2[j]*age2[i,j,k] +
        beta3[j]*age3[i,j,k] + beta4[j]*age4[i,j,k] + beta5[j]*age5[i,j,k] +
         gamma1[k]*male[i,j,k] + gamma2[k]*female[i,j,k] + b[i]
      }
    }
  }
  
  # priors for fixed effects
  # intercepts
  alpha ~ dnorm(mu.int, tau.int)
  
  for (j in 1:n_age) {
    beta1[j] ~ dnorm(mu.beta, tau.beta) # age effects
    beta2[j] ~ dnorm(mu.beta, tau.beta)
    beta3[j] ~ dnorm(mu.beta, tau.beta)
    beta4[j] ~ dnorm(mu.beta, tau.beta)
    beta5[j] ~ dnorm(mu.beta, tau.beta)
  }
  
  for (k in 1:n_sex) {
    gamma1[k] ~ dnorm(mu.gamma, tau.gamma) # gender effects
    gamma2[k] ~ dnorm(mu.gamma, tau.gamma) # gender effects
  }
  
  mu.int ~ dnorm(0, 0.001) # Hyperparameter for random intercepts
  tau.int <- 1 / (sigma.int * sigma.int)
  sigma.int ~ dunif(0, 10)
  
  mu.beta ~ dnorm(0, 0.001)  # Hyperparameter for random slopes
  tau.beta <- 1 / (sigma.beta * sigma.beta)
  sigma.beta ~ dunif(0, 10)
  
  mu.gamma ~ dnorm(0, 0.001)  # Hyperparameter for random slopes
  #mu.gamma2 ~ dnorm(0, 0.001)
  tau.gamma <- 1 / (sigma.gamma * sigma.gamma)
  #tau.gamma2 <- 1 / (sigma.gamma2 * sigma.gamma2)
  sigma.gamma ~ dunif(0, 10)
  #sigma.gamma2 ~ dunif(0, 10)
  
  # Random effects for municipalities
  for (i in 1:n_munic) {
    b[i] ~ dnorm(0, 1.0E-6)
  }
  
  # Prior for precision and standard deviation
  #sigma ~ dunif(0,10)
	#tau <- 1 / (sigma*sigma)
  
}
", fill = TRUE)

sink()


# Initial values
inits2 <- function() {
  list(
    alpha = 0,
    beta1 = rep(0, n_age),
    beta2 = rep(0, n_age),
    beta3 = rep(0, n_age),
    beta4 = rep(0, n_age),
    beta5 = rep(0, n_age),
    gamma1 = rep(0, n_sex),
    gamma2 = rep(0, n_sex),
    b = rep(0, n_munic),
    mu.int = rnorm(1, 0, 1),
    mu.beta = rnorm(1, 0, 1),
    mu.gamma = rnorm(1, 0, 1)
  )
}

params <- c("pi", "alpha", "beta1", "beta2", "beta3", "beta4", "beta5", "gamma1", 
            "gamma2", "mu.int", "sigma.int", "mu.beta",
            "sigma.beta", "mu.gamma", "sigma.gamma")


# Run the model
mod2.fit <- jags(
  data = datajags2,
  inits = inits2,
  parameters.to.save = params,
  model.file = "model2.txt",
  n.chains = 3,
  n.iter = 100000,
  n.burnin=20000,
  jags.seed = 123,
  quiet = FALSE
)

bayes2mcmc <- as.mcmc(mod2.fit)
summary(bayes2mcmc)

autocorr.plot(bayes2mcmc)
geweke.diag(bayes2mcmc)
heidel.diag(bayes2mcmc)

caterplot(bayes2mcmc, parms = c("pi"))

