# load libraries ------------------------------------
library(tidyverse)
library(R2jags)
library(coda)
library(readxl)
library(pacman)
library(MCMCvis)
library(mcmcplots)

# Leo's draft ---------------------------------------

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
  n_munic_arr = array(1:n_munic, dim = c(300, 5, 2)),
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

## winbugs code ----

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
       
       logit(pi[i,j,k]) <- beta1 + beta2*age1[i,j,k] + beta3*age2[i,j,k] + 
        beta4*age3[i,j,k] + beta5*age4[i,j,k] + beta6*sex_num[i,j,k] + 
        b[n_munic_arr[i, j, k]]
       
        
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
  
  # Random effects for municipalities
  for (i in 1:n_munic) {
    b[i] ~ dnorm(0.0,tau.b)
  }
  # priors for random effects
  sigma.b ~ dunif(0,100)
	sigma.2b  <- pow(sigma.b, 2)
	tau.b <- pow(sigma.2b, -1)
	
	# priors for betas
	beta1 ~ dnorm(0, 0.001)
	beta2 ~ dnorm(0, 0.001)
	beta3 ~ dnorm(0, 0.001)
	beta4 ~ dnorm(0, 0.001)
	beta5 ~ dnorm(0, 0.001)
	beta6 ~ dnorm(0, 0.001)
	
}
", fill = TRUE)

sink()

inits2 <- list(
  list(
    beta1 = 0,
    beta2 = 0,
    beta3 = 0,
    beta4 = 0,
    beta5 = 0,
    beta6 = 0,
    sigma.b = 1
  ), list(
    beta1 = 0,
    beta2 = 0,
    beta3 = 0,
    beta4 = 0,
    beta5 = 0,
    beta6 = 0,
    sigma.b = 10
  ), list(
    beta1 = 0,
    beta2 = 0,
    beta3 = 0,
    beta4 = 0,
    beta5 = 0,
    beta6 = 0,
    sigma.b = 5
  )
)

# update initial values for better convergence

params <- c("alpha", "beta1", "beta2", "sigma.b") #, 'B', "mu", "gamma1",
#"mu.int", "sigma.int", "mu.beta", "sigma.beta", 'b'

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
traceplot(mcmc_samples[, c("alpha", "beta1[1]", "beta1[2]", "beta1[3]", 
                           "beta1[4]", "beta1[5]", "gamma1[1]", "gamma1[2]", 
                           "sigma.b", "mu.gamma", "b[2]")])

autocorr.plot(bayes2mcmc)
geweke.diag(bayes2mcmc)
heidel.diag(bayes2mcmc)

caterplot(bayes2mcmc, parms = c("pi"))

# the version works ----
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

## winbugs code ----

sink("model2.txt")
cat("
# ...
model {
  for (i in 1:n_munic) {        # Loop over municipalities
    for (j in 1:n_age) {        # Loop over age groups
      for (k in 1:n_sex) {      # Loop over genders

        Yiag[i,j,k] ~ dbin(pi[i,j,k], Niag[i,j,k])  # Likelihood
        
        mu[i,j,k] <- alpha + beta1[j]*age1[i,j,k] + beta1[j]*age2[i,j,k] +
        beta1[j]*age3[i,j,k] + beta1[j]*age4[i,j,k] + beta1[j]*age5[i,j,k] +
         gamma1[k]*male[i,j,k]
         
        logit(pi[i,j,k]) <- mu[i,j,k] + b[i]
      }
    }
  }
  
  # priors for fixed effects
  # intercepts
  alpha ~ dnorm(0.0,1.0E-6)
  
  for (j in 1:n_age) {
    beta1[j] ~ dnorm(mu.beta, tau.beta) # age effects
  }
  
  for (k in 1:n_sex) {
    gamma1[k] ~ dnorm(mu.gamma,tau.gamma) # gender effects
  }
  
  mu.beta ~ dnorm(0, 0.001)  # Hyperparameter for random slopes
  tau.beta <- 1 / (sigma.beta * sigma.beta)
  sigma.beta ~ dunif(0, 10)
  
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
    sigma.b = runif(1, 0, 10), 
    mu.beta = rnorm(1, 0, 1),
    sigma.beta = runif(1, 0, 10),
    mu.gamma = rnorm(1, 0, 1),
    sigma.gamma = runif(1, 0, 10)
  )
}

params <- c("pi", "alpha", "beta1", "gamma1",
            'b', "mu.beta", "sigma.beta",
            "mu.gamma", "sigma.gamma", "sigma.b")

# Run the model
mod2.fit <- jags(
  data = datajags2,
  inits = inits2,
  parameters.to.save = params,
  model.file = "model2.txt",
  n.chains = 3,
  n.iter = 2000,
  n.burnin = 100,
  jags.seed = 123,
  quiet = FALSE
)

## trace plot ---
mcmc_samples <- as.mcmc(mod2.fit)
traceplot(mcmc_samples[, c("alpha", "beta1[1]", "beta1[2]", "beta1[3]", 
                           "beta1[4]", "beta1[5]", "gamma1[1]", "gamma1[2]", 
                           "sigma.b", "mu.gamma", "b[2]")])

# Edd's book ver ----
# this works too!
DDKdata <- read.csv("data/bda_data.csv")

DDKdata <- DDKdata  |> 
  arrange(sex, age, naam)

datajags2 <- list(
  n_munic = length(unique(DDKdata$naam)),
  n_munic_arr = array(1:n_munic, dim = c(300, 5, 2)),
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

## winbugs code ----

sink("model2.txt")
cat("
model {
# Priors
  for (i in 1:n_munic){
    alpha[i] ~ dnorm(mu.int, tau.int) # Intercepts
}
  mu.int ~ dnorm(0, 0.001) # Hyperparameter for random intercepts
  tau.int <- 1 / (sigma.int * sigma.int)
  sigma.int ~ dunif(0, 10)
  
  for (j in 1:n_age) {
    beta1[j] ~ dnorm(0, 0.001) # age effects
  }
  
  for (k in 1:n_sex) {
    gamma1[k] ~ dnorm(0, 0.001) # sex effects
  }
  
  # Binomial likelihood
  for (i in 1:n_munic) {        # Loop over municipalities
    for (j in 1:n_age) {        # Loop over age groups
      for (k in 1:n_sex) {      # Loop over genders
  Yiag[i,j,k] ~ dbin(pi[i,j,k], Niag[i,j,k])
  logit(pi[i,j,k]) <- alpha[n_munic_arr[i, j, k]] + beta1[j]*age1[i,j,k] + 
  beta1[j]*age2[i,j,k] + beta1[j]*age3[i,j,k] + beta1[j]*age4[i,j,k] + 
  beta1[j]*age5[i,j,k] + gamma1[k]*male[i,j,k] + gamma1[k]*female[i,j,k]
      }
    }
  }
}
", fill = TRUE)
sink()

# Initial values
inits2 <- function() {
  list(
    alpha = rnorm(n_munic, 0, 2),
    beta1 = rnorm(n_age, 1, 1),
    gamma1 = rnorm(n_sex, 1, 1),
    mu.int = rnorm(1, 0, 1)
  )
}

params <- c("alpha", "beta1", "gamma1", "mu.int", "sigma.int")

# Run the model
mod2.fit <- jags(
  data = datajags2,
  inits = inits2,
  parameters.to.save = params,
  model.file = "model2.txt",
  n.chains = 3,
  n.iter = 2000,
  n.burnin = 100,
  jags.seed = 123,
  quiet = FALSE
)

## trace plot ---
mcmc_samples <- as.mcmc(mod2.fit)
traceplot(mcmc_samples[, c("alpha[1]", "beta1[1]", "beta1[2]", "beta1[3]", 
                           "beta1[4]", "beta1[5]", "gamma1[1]", "gamma1[2]")])

# Edd's book modified ----
# this doesn't work

DDKdata <- read.csv("data/bda_data.csv")

DDKdata <- DDKdata  |> 
  arrange(sex, age, naam)

datajags2 <- list(
  n_munic = length(unique(DDKdata$naam)),
  n_munic_arr = array(1:n_munic, dim = c(300, 5, 2)),
  n_age = length(unique(DDKdata$age)),
  n_sex = length(unique(DDKdata$sex)),
  Niag = array(c(DDKdata$invited), dim = c(300, 5, 2)), 
  Yiag = array(c(DDKdata$participant), dim = c(300, 5, 2)),
  male = array(c(DDKdata$male), dim = c(300, 5, 2)),
  age1 = array(c(DDKdata$age1), dim = c(300, 5, 2)),
  age2 = array(c(DDKdata$age2), dim = c(300, 5, 2)),
  age3 = array(c(DDKdata$age3), dim = c(300, 5, 2)),
  age4 = array(c(DDKdata$age4), dim = c(300, 5, 2))
)

## winbugs code ----

sink("model2.txt")
cat("
model {
# Priors
  for (i in 1:n_munic){
    alpha[i] ~ dnorm(mu.int, tau.int) # Intercepts
}
  mu.int ~ dnorm(0, 0.001) # Hyperparameter for random intercepts
  tau.int <- 1 / (sigma.int * sigma.int)
  sigma.int ~ dunif(0, 10)
  
  beta ~ dnorm(0, 0.001) # age effects
  gamma1 ~ dnorm(0, 0.001) # sex effects
  
  # Binomial likelihood
  for (i in 1:n_munic) {        # Loop over municipalities
    for (j in 1:n_age) {        # Loop over age groups
      for (k in 1:n_sex) {      # Loop over genders
  Yiag[i,j,k] ~ dbin(pi[i,j,k], Niag[i,j,k])
  logit(pi[i,j,k]) <- alpha[n_munic_arr[i, j, k]] + beta*age1[i,j,k] + 
  beta*age2[i,j,k] + beta*age3[i,j,k] + beta*age4[i,j,k] + gamma1*male[i,j,k]
      }
    }
  }
}
", fill = TRUE)
sink()

# Initial values
inits2 <- function() {
  list(
    alpha = rnorm(n_munic, 0, 2),
    beta = rnorm(1, 1, 1),
    gamma1 = rnorm(1, 1, 1),
    mu.int = rnorm(1, 0, 1)
  )
}

params <- c("alpha", "beta", "gamma1", "mu.int", "sigma.int")

# Run the model
mod2.fit <- jags(
  data = datajags2,
  inits = inits2,
  parameters.to.save = params,
  model.file = "model2.txt",
  n.chains = 3,
  n.iter = 2000,
  n.burnin = 100,
  jags.seed = 123,
  quiet = FALSE
)

## trace plot ---
mcmc_samples <- as.mcmc(mod2.fit)
traceplot(mcmc_samples[, c("alpha[1]", "beta1[1]", "beta1[2]", "beta1[3]", 
                           "beta1[4]", "beta1[5]", "gamma1[1]", "gamma1[2]")])

# edd's modified 2 ----
DDKdata <- read.csv("data/bda_data.csv")

DDKdata <- DDKdata  |> 
  arrange(sex, age, naam)

datajags2 <- list(
  n_munic = length(unique(DDKdata$naam)),
  n_munic_arr = array(1:n_munic, dim = c(300, 5, 2)),
  n_age = length(unique(DDKdata$age)),
  n_sex = length(unique(DDKdata$sex)),
  Niag = array(c(DDKdata$invited), dim = c(300, 5, 2)), 
  Yiag = array(c(DDKdata$participant), dim = c(300, 5, 2)),
  male = array(c(DDKdata$male), dim = c(300, 5, 2)),
  age1 = array(c(DDKdata$age1), dim = c(300, 5, 2)),
  age2 = array(c(DDKdata$age2), dim = c(300, 5, 2)),
  age3 = array(c(DDKdata$age3), dim = c(300, 5, 2)),
  age4 = array(c(DDKdata$age4), dim = c(300, 5, 2))
)

## winbugs code ----

sink("model2.txt")
cat("
model {
# Priors
  for (i in 1:n_munic){
    alpha[i] ~ dnorm(mu.int, tau.int) # Intercepts
}
  mu.int ~ dnorm(0, 0.001) # Hyperparameter for random intercepts
  tau.int <- 1 / (sigma.int * sigma.int)
  sigma.int ~ dunif(0, 10)
  
  for (j in 1:n_age) {
    beta1[j] ~ dnorm(0, 0.001) # age effects
  }
  
  for (k in 1:n_sex) {
    gamma1[k] ~ dnorm(0, 0.001) # sex effects
  }
  
  # Binomial likelihood
  for (i in 1:n_munic) {        # Loop over municipalities
    for (j in 1:n_age) {        # Loop over age groups
      for (k in 1:n_sex) {      # Loop over genders
  Yiag[i,j,k] ~ dbin(pi[i,j,k], Niag[i,j,k])
  logit(pi[i,j,k]) <- alpha[n_munic_arr[i, j, k]] + beta1[j]*age1[i,j,k] + 
  beta1[j]*age2[i,j,k] + beta1[j]*age3[i,j,k] + beta1[j]*age4[i,j,k] + 
  gamma1[k]*male[i,j,k]
      }
    }
  }
}
", fill = TRUE)
sink()

# Initial values
inits2 <- function() {
  list(
    alpha = rnorm(n_munic, 0, 2),
    beta1 = rnorm(n_age, 1, 1),
    gamma1 = rnorm(n_sex, 1, 1),
    mu.int = rnorm(1, 0, 1)
  )
}

params <- c("alpha", "beta1", "gamma1", "mu.int", "sigma.int")

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

## trace plot ---
mcmc_samples <- as.mcmc(mod2.fit)
traceplot(mcmc_samples[, c("alpha[1]", "beta1[1]", "beta1[2]", "beta1[3]", 
                           "beta1[4]", "gamma1[2]", "mu.int", 
                           "sigma.int", "alpha[2]", "alpha[3]", "alpha[4]", 
                           "alpha[5]")])

# Leo's final ---------------------------------------
## data preparation ----
DDKdata <- read.csv("data/bda_data.csv")

DDKdata <- DDKdata  |> 
  arrange(sex, age, naam)

datajags2 <- list(
  n_munic = length(unique(DDKdata$naam)),
  n_munic_arr = array(1:n_munic, dim = c(300, 5, 2)),
  n_age = length(unique(DDKdata$age)),
  n_sex = length(unique(DDKdata$sex)),
  Niag = array(c(DDKdata$invited), dim = c(300, 5, 2)), 
  Yiag = array(c(DDKdata$participant), dim = c(300, 5, 2)),
  male = array(c(DDKdata$male), dim = c(300, 5, 2)),
  age1 = array(c(DDKdata$age1), dim = c(300, 5, 2)),
  age2 = array(c(DDKdata$age2), dim = c(300, 5, 2)),
  age3 = array(c(DDKdata$age3), dim = c(300, 5, 2)),
  age4 = array(c(DDKdata$age4), dim = c(300, 5, 2))
)

## winbugs code ----

sink("model2.txt")
cat("
model {
# Priors
  for (i in 1:n_munic){
    alpha[i] ~ dnorm(mu.int, tau.int) # Intercepts
}
  mu.int ~ dnorm(0, 0.001) # Hyperparameter for random intercepts
  tau.int <- 1 / (sigma.int * sigma.int)
  sigma.int ~ dunif(0, 10)
  
  for (j in 1:n_age) {
    beta1[j] ~ dnorm(0, 0.001) # age effects
  }
  
  for (k in 1:n_sex) {
    gamma1[k] ~ dnorm(0, 0.001) # sex effects
  }
  
  # Binomial likelihood
  for (i in 1:n_munic) {        # Loop over municipalities
    for (j in 1:n_age) {        # Loop over age groups
      for (k in 1:n_sex) {      # Loop over genders
  Yiag[i,j,k] ~ dbin(pi[i,j,k], Niag[i,j,k])
  logit(pi[i,j,k]) <- alpha[n_munic_arr[i, j, k]] + beta1[j]*age1[i,j,k] + 
  beta1[j]*age2[i,j,k] + beta1[j]*age3[i,j,k] + beta1[j]*age4[i,j,k] + 
  gamma1[k]*male[i,j,k]
      }
    }
  }
}
", fill = TRUE)
sink()

# Initial values
inits2 <- function() {
  list(
    alpha = rnorm(n_munic, 0, 2),
    beta1 = rnorm(n_age, 1, 1),
    gamma1 = rnorm(n_sex, 1, 1),
    mu.int = rnorm(1, 0, 1)
  )
}

params <- c("alpha", "beta1", "gamma1", "mu.int", "sigma.int", 'pi')

## Run the model ----
mod2.fit <- jags(
  data = datajags2,
  inits = inits2,
  parameters.to.save = params,
  model.file = "model2.txt",
  n.chains = 3,
  n.iter = 10000,
  n.burnin = 5000,
  jags.seed = 123,
  quiet = FALSE
)

# print pD and DIC ---------------------------------------
options(max.print=999999)
sink("model2res1.txt")
print(mod2.fit)
sink()

options(max.print=999999)
sink("model2res1.txt")
bayes2mcmc <- as.mcmc(mod2.fit)
summary(bayes2mcmc)
sink()

summary(bayes2mcmc)
## diagnostics ----------------------------------------------
mcmc_samples <- as.mcmc(mod2.fit)

png("pictures/mod2trace.png", width = 18, 
    height = 15, units = "cm", res = 300)
#traceplot(mcmc_samples[, c("alpha[1]", "beta1[1]", "beta1[2]", "beta1[3]", 
#                           "beta1[4]", "gamma1[2]", "mu.int", 
#                           "sigma.int", "alpha[2]", "alpha[3]", "alpha[4]", 
#                           "alpha[5]")])
dev.off()

## library coda -----------------------------------------------
library(coda)

## autocorrelation
autocorr.diag(as.mcmc(mod2.fit))
autocorr.plot(as.mcmc(mod2.fit))

png("pictures/mod2rmean_sigmaint.png", width = 18, 
    height = 10, units = "cm", res = 300)
rmeanplot(as.mcmc(mod2.fit), parms = "sigma.int")
dev.off()

## geweke diag - works
png("pictures/mod2geweke_sigmaint.png", width = 18, 
    height = 10, units = "cm", res = 300)
geweke.plot(as.mcmc(mod2.fit), parms = "sigma.int")
dev.off()

# Heidel diag
heidel.diag(as.mcmc(mod2.fit))

# BGR diagnostics
gelman.diag(as.mcmc(mod2.fit))
gelman.plot(as.mcmc(mod2.fit))

# Raftery-Lewis
raftery.diag(as.mcmc(mod2.fit))

# HW diagnostic
sink("model2heidel.txt")
heidel.diag(as.mcmc(mod2.fit))
sink()

### library ggmcmc -----------------------------------------------
library(ggmcmc)
bayes2mcmcggs <- ggs(bayes2mcmc)

png("pictures/mod2trace_muint.png", width = 18, 
    height = 10, units = "cm", res = 300)
ggs_traceplot(bayes2mcmcggs, family = "mu.int")
dev.off()

png("pictures/mod2autocorr_gamma.png", width = 18, 
    height = 10, units = "cm", res = 300)
ggs_autocorrelation(bayes2mcmcggs, family = "gamma1")
dev.off()

png("pictures/mod2geweke_sigmaint.png", width = 18, 
    height = 10, units = "cm", res = 300)
ggs_geweke(bayes2mcmcggs, family = "sigma.int")
dev.off()

png("pictures/mod2gbr_gamma.png", width = 18, 
    height = 10, units = "cm", res = 300)
ggs_grb(bayes2mcmcggs, family = "gamma1")
dev.off()


print(ggs_diagnostics(bayes2mcmcggs, family = "pi"), n=500)

#------------------------------------------------------------------
### Extract posterior samples for pi

# extract chains
ex2 <- MCMCchains(bayes2mcmc, params = 'pi')

# Compute P(pi_i < 0.30) for each column
(probs2 <- apply(ex2, 2, function(x) mean(x < 0.30)))

# Find columns where P(pi_i < 0.30) > 0.9
select2 <- which(probs2 > 0.9)

# Display results
list(
  Columns = select2,
  Probabilities = probs[select2]
)


### --------------------------------------------------------
library(BayesFactor)
load(file='mod.fit')

compare(mod.fit, mod2.fit) # doesnt work


####################### outliers detection -------------------------------------------------------
# Load posterior samples from JAGS
library(coda)

jags2 <- jags.model(file="model2.txt",
                   data = datajags2,
                   inits = inits2,
                   n.chains = 3)

update(jags2, 5000) # burn-in period

samples <- coda.samples(model = jags2,
                          variable.names = params,
                          n.iter=10000, 
                          thin=5)


### ############## HIERARCHICAL CENTERING ----------------

sink("model22.txt")
cat("
model {
  # Priors
    for (i in 1:n_munic){
      alpha[i] ~ dnorm(mu.int, tau.int) # Intercepts
  }
    mu.int ~ dnorm(0, 0.001) # Hyperparameter for random intercepts
    tau.int <- 1 / (sigma.int * sigma.int)
    sigma.int ~ dunif(0, 10)
  
  for (j in 1:n_age) {
    beta1[j] ~ dnorm(0, 0.001) # age effects
  }
  
  for (k in 1:n_sex) {
    gamma1[k] ~ dnorm(0, 0.001) # sex effects
  }
  
  # Binomial likelihood
  for (i in 1:n_munic) {        # Loop over municipalities
    for (j in 1:n_age) {        # Loop over age groups
      for (k in 1:n_sex) {      # Loop over genders
      
        Yiag[i,j,k] ~ dbin(pi[i,j,k], Niag[i,j,k])
        
        m[i, j, k] <- alpha[n_munic_arr[i, j, k]] + beta1[j]*age1[i,j,k] + 
        beta1[j]*age2[i,j,k] + beta1[j]*age3[i,j,k] + beta1[j]*age4[i,j,k] + 
        gamma1[k]*male[i,j,k]
      
        n[i, j, k] <- m[i, j, k] + alpha[n_munic_arr[i, j, k]]
        
        logit(pi[i,j,k]) <- n[i, j, k]
      
        #n[i, j, k] ~ dnorm(m[i, j, k], tau)
      }
    }
  }
}
", fill = TRUE)
sink()

inits2 <- function() {
  list(
    alpha = rnorm(n_munic, 0, 2),
    beta1 = rnorm(n_age, 1, 1),
    gamma1 = rnorm(n_sex, 1, 1),
    mu.int = rnorm(1, 0, 1)
  )
}

params <- c("alpha", "beta1", "gamma1", "mu.int", "sigma.int")

## Run the model ----
mod22.fit <- jags(
  data = datajags2,
  inits = inits2,
  parameters.to.save = params,
  model.file = "model22.txt",
  n.chains = 3,
  n.iter = 10000,
  n.burnin = 5000,
  jags.seed = 123,
  quiet = FALSE
)

options(max.print=999999)
sink("model22info.txt")
print(mod22.fit)
sink()

options(max.print=999999)
sink("model22res.txt")
bayes22mcmc <- as.mcmc(mod22.fit)
summary(bayes22mcmc)
sink()

#----------------------------------------------

