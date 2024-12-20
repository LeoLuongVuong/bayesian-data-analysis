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
  #gender = array(c(DDKdata$sex), dim = c(300, 5, 2)),
  agec = array(c(DDKdata$agec), dim = c(300, 5, 2))
)

n_munic = length(unique(DDKdata$naam))
n_age = length(unique(DDKdata$age))
n_sex = length(unique(DDKdata$sex))
agec = array(c(DDKdata$agec), dim = c(300, 5, 2))
gender = array(c(DDKdata$sex), dim = c(300, 5, 2))

dim(agec)

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
        logit(pi[i,j,k]) <- alpha + beta[j]*agec[i,j,k] + b[i]
      }
    }
  }
  
  # priors for fixed effects
  # intercepts
  alpha ~ dnorm(mu.int, tau.int)
  
  for (j in 1:n_age) {
    beta[j] ~ dnorm(mu.beta, tau.beta) # age effects
  }
  
  for (k in 1:n_sex) {
    gamma[k] ~ dnorm(0, 1.0E-6) # gender effects
  }
  
  mu.int ~ dnorm(0, 0.001) # Hyperparameter for random intercepts
  tau.int <- 1 / (sigma.int * sigma.int)
  sigma.int ~ dunif(0, 10)
  
  mu.beta ~ dnorm(0, 0.001)  # Hyperparameter for random slopes
  tau.beta <- 1 / (sigma.beta * sigma.beta)
  sigma.beta ~ dunif(0, 10)
  
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
    beta = rep(0, n_age),
    gamma = rep(0, n_sex),
    b = rep(0, n_munic),
    mu.int = rnorm(1, 0, 1),
    mu.beta = rnorm(1, 0, 1)
  )
}

params <- c("pi","alpha", "beta", "mu.int", "sigma.int", "mu.beta", "sigma.beta")


# Run the model
mod2.fit <- jags(
  data = datajags2,
  inits = inits2,
  parameters.to.save = params,
  model.file = "model2.txt",
  n.chains = 3,
  n.iter = 2000,
  n.burnin=1000,
  jags.seed = 123,
  quiet = FALSE
)

bayes2mcmc <- as.mcmc(mod2.fit)
summary(bayes2mcmc)

caterplot(bayes2mcmc, parms = c("pi"))

