# load libraries ------------------------------------
library(tidyverse)
library(R2jags)
library(coda)
library(readxl)
library(pacman)
library(MCMCvis)
library(mcmcplots)

# data manipulation --------------------------------

bda_data <- read.csv("data/bda_data.csv")

bda_data <- bda_data %>%
  arrange(sex, age, naam)

# sex and age to factors
bda_data$sex <- as.factor(bda_data$sex)
bda_data$age <- as.factor(bda_data$age)
bda_data$age_num <- as.numeric(bda_data$age)
bda_data$sex_num <- ifelse(bda_data$sex == "female", 0, 1)

datajags2 <- list(
  n_munic = length(unique(bda_data$naam)),
  n_age = length(unique(bda_data$age)),
  n_sex = length(unique(bda_data$sex)),
  Niag = array(c(bda_data$invited), dim = c(300, 5, 2)), 
  Yiag = array(c(bda_data$participant), dim = c(300, 5, 2)),
  sex_num = array(c(bda_data$sex_num), dim = c(300, 5, 2)),
  age1 = array(c(bda_data$age1), dim = c(300, 5, 2)),
  age2 = array(c(bda_data$age2), dim = c(300, 5, 2)),
  age3 = array(c(bda_data$age3), dim = c(300, 5, 2)),
  age4 = array(c(bda_data$age4), dim = c(300, 5, 2))
  
)

n_munic = length(unique(bda_data$naam))
n_age = length(unique(bda_data$age))
n_sex = length(unique(bda_data$sex))
# agec = array(c(bda_data$agec), dim = c(300, 5, 2))
# gender = array(c(bda_data$sex), dim = c(300, 5, 2))

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
       
       logit(pi[i,j,k]) <- beta[1] + beta[2]*age1[i,j,k] + beta[3]*age2[i,j,k] + 
        
        beta[4]*age3[i,j,k] + beta[5]*age4[i,j,k] + beta[6]*sex_num[i,j,k] + b[i]
       
        
        # link function:logit of participation rate
        
        # mu[i, j, k] <- beta[1] + beta[2]*age1[i,j,k] + beta[3]*age2[i,j,k] + 
        # beta[4]*age3[i,j,k] + beta[5]*age4[i,j,k] + beta[6]*sex_num[i,j,k]
         
        # logit(pi[i,j,k]) <- mu[i, j, k] + b[i]
         
         # remove gamma2[k]*female[i,j,k] since male is a dummy in itself
         # give all age the same beta
         
         # compare predicted vs observed prob to examine shrinkage
          # B[i, j, k] <- (exp(mu[i,j,k]) / (1 + exp(mu[i,j,k]))) / (Yiag[i,j,k]/Niag[i,j,k])
      }
    }
  }
  
  # priors for fixed effects
  # intercepts
  #alpha ~ dnorm(0.0,1.0E-2)
  
  #gamma1 ~ dnorm(0.0,1.0E-2)
  
  # Random effects for municipalities
  for (i in 1:n_munic) {
    b[i] ~ dnorm(0.0,tau.b)
  }
  # priors for random effects
  sigma.b ~ dunif(0,100)
	sigma.2b  <- pow(sigma.b, 2)
	tau.b <- pow(sigma.2b, -1)
	
	for (l in 1:6) {
     # age effects
   beta[l]  ~ dnorm(0.0,1.0E-4)
  }
}
", fill = TRUE)

sink()
  
inits2 <- list(
  list(
    #alpha = rnorm(1, 0, 1),
    beta = rnorm(6, 0, 1),
    #gamma1 = rnorm(1, 0, 1),
    #b = rnorm(n_munic, 0, 1),
    sigma.b = runif(1, 0, 10)
  ), list(
    #alpha = rnorm(1, 0, 1),
    beta = rnorm(6, 0, 1),
    #gamma1 = rnorm(1, 0, 1),
    #b = rnorm(n_munic, 0, 1),
    sigma.b = runif(1, 0, 10)
  ), list(
    #alpha = rnorm(1, 0, 1),
    beta = rnorm(6, 0, 1),
    #gamma1 = rnorm(1, 0, 1),
    #b = rnorm(n_munic, 0, 1),
    sigma.b = runif(1, 0, 10)
  )
)

# update initial values for better convergence

params <- c("pi", "alpha", "beta", 
           'b', "sigma.b") #, 'B', "mu", "gamma1",
#"mu.int", "sigma.int", "mu.beta", "sigma.beta",

# Run the model
mod2.fit <- jags(
  data = datajags2,
  inits = inits2,
  parameters.to.save = params,
  model.file = "model2.txt",
  n.chains = 3,
  n.iter = 5000,
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

