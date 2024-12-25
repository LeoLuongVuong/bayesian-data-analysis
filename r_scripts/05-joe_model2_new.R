library(tidyverse)
library(R2jags)
library(coda)
library(readxl)
library(pacman)
library(MCMCvis)
library(mcmcplots)

DDKdata <- read.csv("data/bda_data.csv")

DDKdata <- DDKdata  |> 
  arrange(sex, age, naam)

datajags2 <- list(
  n_munic = length(unique(DDKdata$naam)),
  n_age = length(unique(DDKdata$age)),
  n_sex = length(unique(DDKdata$sex)),
  Niag = array(c(DDKdata$invited), dim = c(300, 5, 2)), 
  Yiag = array(c(DDKdata$participant), dim = c(300, 5, 2)),
  male = array(c(DDKdata$male), dim = c(300, 5, 2)),
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

# Q2 - Leo ----

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
        
        mu[i,j,k] <- alpha + beta1[j]*age1[i,j,k] + beta1[j]*age2[i,j,k] +
        beta1[j]*age3[i,j,k] + beta1[j]*age4[i,j,k] + beta1[j]*age5[i,j,k] +
         gamma1[k]*male[i,j,k]
         
        logit(pi[i,j,k]) <- mu[i,j,k] + b[i]
         
         # remove gamma2[k]*female[i,j,k] since male is a dummy in itself
         # give all age the same beta
         
         # compare predicted vs observed prob to examine shrinkage
          #B[i, j, k] <- (exp(mu[i,j,k]) / (1 + exp(mu[i,j,k]))) / (Yiag[i,j,k]/Niag[i,j,k])
      }
    }
  }
  
  # priors for fixed effects
  # intercepts
  alpha ~ dnorm(0.0,1.0E-6)
  
  for (j in 1:n_age) {
    beta1[j] ~ dnorm(mu.beta,tau.beta) # age effects
    # remove beta 2 to beta 5 j
  }
  
  for (k in 1:n_sex) {
    gamma1[k] ~ dnorm(mu.gamma,tau.gamma) # gender effects
    # remove gamma2[k]
  }
  
  # mu.int ~ dnorm(0, 0.001) # Hyperparameter for random intercepts
  # tau.int <- 1 / (sigma.int * sigma.int)
  # sigma.int ~ dunif(0, 10)
  
  mu.beta ~ dnorm(0, 0.001)  # Hyperparameter for random slopes
  tau.beta <- 1 / (sigma.beta * sigma.beta)
  sigma.beta ~ dunif(0, 10)
  # 
  mu.gamma ~ dnorm(0, 0.001)  # Hyperparameter for random slopes
  tau.gamma <- 1 / (sigma.gamma * sigma.gamma)
  sigma.gamma ~ dunif(0, 10)
  
  # Random effects for municipalities
  for (i in 1:n_munic) {
    b[i] ~ dnorm(0.0,tau.b)
  }
  tau.b <- 1 / (sigma.b * sigma.b)
  sigma.b ~ dunif(0, 10)
  
  
  
}
", fill = TRUE)

sink()


# Initial values
inits2 <- function() {
  list(
    alpha = 0,
    beta1 = rep(0, n_age),
    gamma1 = rep(0, n_sex),
    b = rep(0, n_munic),
    #mu.int = rnorm(1, 0, 1),
    sigma.b = runif(1, 0, 10), #, n_munic
    mu.beta = rnorm(1, 0, 1),
    sigma.beta = runif(1, 0, 10),
    mu.gamma = rnorm(1, 0, 1),
    sigma.gamma = runif(1, 0, 10)
  )
}

params <- c("pi", "alpha", "beta1", "gamma1",
           'b', "mu.beta", "sigma.beta",
           "mu.gamma", "sigma.gamma", "sigma.b") #, 'B', "mu", 
#"mu.int", "sigma.int", "mu.beta", "sigma.beta",


# Run the model
mod2.fit <- jags(
  data = datajags2,
  inits = inits2,
  parameters.to.save = params,
  model.file = "model2.txt",
  n.chains = 3,
  n.iter = 10000,
  n.burnin = 100,
  jags.seed = 123,
  quiet = FALSE
)

## trace plot ---
mcmc_samples <- as.mcmc(mod2.fit)
traceplot(mcmc_samples[, c("alpha", "beta1[1]", "gamma1[1]", "sigma.b")])

autocorr.plot(bayes2mcmc)
geweke.diag(bayes2mcmc)
heidel.diag(bayes2mcmc)

caterplot(bayes2mcmc, parms = c("pi"))

