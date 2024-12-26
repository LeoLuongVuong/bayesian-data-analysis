# load libraries ------------------------------------
library(tidyverse)
library(R2jags)
library(coda)
library(readxl)
library(pacman)
library(MCMCvis)
library(mcmcplots)

# bda_data_long <- read.csv("data/bda_data.csv")
# 
# bda_data_long <- bda_data_long  |> 
#   arrange(sex, age, naam)

# data manipulation --------------------------------

bda_data_long <- read.csv("data/bda_data_long.csv")

bda_data_long <- bda_data_long %>%
  arrange(sex, age, naam)

# sex and age to factors
bda_data_long$sex <- as.factor(bda_data_long$sex)
bda_data_long$age <- as.factor(bda_data_long$age)
bda_data_long$age_num <- as.numeric(bda_data_long$age)
bda_data_long$sex_num <- ifelse(bda_data_long$sex == "female", 0, 1)

datajags2 <- list(
  n_munic = length(unique(bda_data_long$naam)),
  n_age = length(unique(bda_data_long$age)),
  n_sex = length(unique(bda_data_long$sex)),
  Niag = array(c(bda_data_long$invited), dim = c(300, 5, 2)), 
  Yiag = array(c(bda_data_long$participant), dim = c(300, 5, 2)),
  sex_num = array(c(bda_data_long$sex_num), dim = c(300, 5, 2)),
  age_num = array(c(bda_data_long$age_num), dim = c(300, 5, 2))
)

n_munic = length(unique(bda_data_long$naam))
n_age = length(unique(bda_data_long$age))
n_sex = length(unique(bda_data_long$sex))
# agec = array(c(bda_data_long$agec), dim = c(300, 5, 2))
# gender = array(c(bda_data_long$sex), dim = c(300, 5, 2))

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
        
        mu[i,j,k] <- alpha + beta1*age_num[i,j,k] + gamma1*sex_num[i,j,k]
         
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
  alpha ~ dnorm(0.0,1.0E-2)
  
  beta1 ~ dnorm(0.0,1.0E-2)
  
  gamma1 ~ dnorm(0.0,1.0E-2)
  
  # for (j in 1:n_age) {
  #    # age effects
  #   # remove beta 2 to beta 5 j
  # }
  
  # for (k in 1:n_sex) {
  #    # gender effects
  #   # remove gamma2[k]
  # }
  
  # mu.int ~ dnorm(0, 0.001) # Hyperparameter for random intercepts
  # tau.int <- 1 / (sigma.int * sigma.int)
  # sigma.int ~ dunif(0, 10)
  
  # mu.beta ~ dnorm(0, 0.001)  # Hyperparameter for random slopes
  # tau.beta <- 1 / (sigma.beta * sigma.beta)
  # sigma.beta ~ dunif(0, 10)
  # # 
  # mu.gamma ~ dnorm(0, 0.001)  # Hyperparameter for random slopes
  # tau.gamma <- 1 / (sigma.gamma * sigma.gamma)
  # sigma.gamma ~ dunif(0, 10)
  
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
inits2 <- list(
  list(
    alpha = 0.5,
    beta1 = 0.5,
    gamma1 = 0,
    b = rep(0, n_munic),
    sigma.b = runif(1, 0, 10) #, n_munic
    # mu.beta = rnorm(1, 0, 1),
    # sigma.beta = runif(1, 0, 10),
    # mu.gamma = rnorm(1, 0, 1),
    # sigma.gamma = runif(1, 0, 10)
  ),
  list(
    alpha = -0.5,
    beta1 = -0.5,
    gamma1 = 0.5,
    b = rep(0, n_munic),
    sigma.b = runif(1, 0, 10)),
list(
  alpha = -1,
  beta1 = -1,
  gamma1 = -0.5,
  b = rep(0, n_munic),
  sigma.b = runif(1, 0, 10))
)

# update initial values for better convergence

params <- c("pi", "alpha", "beta1", "gamma1",
           'b', "sigma.b") #, 'B', "mu", 
#"mu.int", "sigma.int", "mu.beta", "sigma.beta",

# Run the model
mod2.fit <- jags(
  data = datajags2,
  inits = inits2,
  parameters.to.save = params,
  model.file = "model2.txt",
  n.chains = 3,
  n.iter = 2000,
  n.burnin = 500,
  jags.seed = 123,
  quiet = FALSE
)

mod3.fit <- jags(
  data = datajags2,
  inits = inits2,
  parameters.to.save = params,
  model.file = "model2.txt",
  n.chains = 3,
  n.iter = 50000,
  n.burnin = 25000,
  jags.seed = 123,
  quiet = FALSE
)

## trace plot ---
mcmc_samples <- as.mcmc(mod2.fit)
traceplot(mcmc_samples[, c("alpha", "beta1[1]", "beta1[2]", "gamma1[1]", 
                           "sigma.b", "mu.gamma", "b[2]")])

autocorr.plot(bayes2mcmc)
geweke.diag(bayes2mcmc)
heidel.diag(bayes2mcmc)

caterplot(bayes2mcmc, parms = c("pi"))

