library(R2jags)
library(coda)
library(readxl)
library(pacman)
library(MCMCvis)

DDKdata <- read.csv("bda_data_long.csv")

datajags2 <- list(
  n_munic = length(unique(DDKdata$naam)),
  n_age = length(unique(DDKdata$age)),
  n_sex = length(unique(DDKdata$sex)),
  Niag = array(c(DDKdata$invited), dim = c(300, 5, 2)),
  Yiag = array(c(DDKdata$participant), dim = c(300, 5, 2))
)

n_munic = length(unique(DDKdata$naam))
n_age = length(unique(DDKdata$age))
n_sex = length(unique(DDKdata$sex))


sink()
sink("model2.txt")

cat("
# ...
model {
  for (i in 1:n_munic) {        # Loop over municipalities
    for (j in 1:n_age) {        # Loop over age groups
      for (k in 1:n_sex) {      # Loop over genders
        # participation rate
        Yiag[i,j,k] ~ dbin(pi[i,j,k], Niag[i,j,k])  # Likelihood
        # logit of participation rate
        logit(pi[i,j,k]) <- alpha + beta[j] + gamma[k] + b[i]
      }
    }
  }
  
  # priors for fixed effects
  alpha ~ dnorm(0.0, 1.0E-6)
  
  for (j in 1:n_age) {
    beta[j] ~ dnorm(0, 1.0E-6) # age effects
  }
  
  for (k in 1:n_sex) {
    gamma[k] ~ dnorm(0, 1.0E-6) # gender effects
  }
  
  # Random effects for municipalities
  for (i in 1:n_munic) {
    b[i] ~ dnorm(0, 1.0E-6)
  }
  
  # Prior for precision and standard deviation
  sigma ~ dunif(0.001,1000)
	tau <- 1 / sqrt(sigma)
  
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
    sigma = 1
  )
}

params <- c("pi", "alpha", "beta", "gamma", "b", "sigma")


# Run the model
mod2.fit <- jags(
  data = datajags,
  inits = inits2,
  parameters.to.save = params,
  model.file = "model2.txt",
  n.chains = 3,
  n.iter = 2000,
  n.burnin=1000,
  jags.seed = 123,
  quiet = FALSE
)


