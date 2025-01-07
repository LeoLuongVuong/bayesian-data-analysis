## load packages ---------------------------
library(rjags)
library(R2jags)
library(coda)
library(readxl)
library(pacman)
library(MCMCvis)
library(mcmcplots)
library(ggmcmc)
library(bayesplot)

pacman::p_load(tidyverse, coda,readxl,janitor)

###-------data preparation  ---------------------------
#### load dataset
bda_data <- read_excel("data/DDK2022.xlsx") %>% 
  janitor::clean_names()

which(is.na(bda_data))

### add participant total
bda_data <-  bda_data %>% 
  mutate(
    totinv = rowSums(select(.,starts_with("invited"))),
    totpart = rowSums(select(.,starts_with("participant"))),
    rate = totpart/totinv
  )
which(is.na(bda_data))

### Reshape data
bda_data_long <- bda_data %>% pivot_longer(
  cols = -c(naam, totinv,totpart, rate),
  names_to = c(".value", "sex", "age"),
  names_pattern = "(invited|participant)_(male|female)_(.*)_yrs"
)

which(is.na(bda_data_long))

## --- write to csv for model 2 
#write.csv(bda_data_long, 'data/bda_data', row.names = FALSE)

bda_data_long <- read.csv("data/bda_data.csv")

### sort bda_data_long in the order of naam, age, sex
bda_data_long <- bda_data_long %>%
  arrange(sex, age, naam)

datajags <- list(
  n_munic = length(unique(bda_data_long$naam)),
  n_age = length(unique(bda_data_long$age)),
  n_sex = length(unique(bda_data_long$sex)),
  Niag = array(c(bda_data_long$invited), dim = c(300, 5, 2)), 
  Yiag = array(c(bda_data_long$participant), dim = c(300, 5, 2))
)

n_munic <- length(unique(bda_data_long$naam))
n_age <- length(unique(bda_data_long$age))
n_sex <- length(unique(bda_data_long$sex))

#### model 1 specification --------------------------------------

sink("model1.txt")
cat("
# ...
model {
  for (i in 1:n_munic) {        # Loop over municipalities
    pi[i] ~ dbeta(1.0, 1.0)         # Prior for participation rate per municipality
    
    for (j in 1:n_age) {        # Loop over age groups
      for (k in 1:n_sex) {      # Loop over genders
        Yiag[i,j,k] ~ dbin(pi[i], Niag[i,j,k])  # Likelihood
      }
    }
  }
}
", fill = TRUE)

sink()

# Initial parameters
my.inits <- list(
  list(pi = rep(0.1, 300)),
  list(pi = rep(0.5, 300)),
  list(pi = rep(0.9, 300))
)


# params to monitor
params <- c("pi")

# Run the model ---- r2jags version ---------------------
mod.fit <- jags(
  data = datajags,
  inits = my.inits, 
  parameters.to.save = c("pi"),
  model.file = "model1.txt",
  n.chains = 3,
  n.iter = 10000,
  n.burnin = 5000,
  jags.seed = 123,
  quiet = FALSE
)

bayes1mcmc <- as.mcmc(mod.fit)

# print pD and DIC ---------------------------------------
options(max.print=999999)
sink("model1results.txt")
print(mod.fit)
sink()

sink("model1res1.txt")
bayes1mcmc <- as.mcmc(mod.fit)
summary(bayes1mcmc)
sink()

##-------convergence checks  --------------------------------------------
traceplot(mod.fit)
rmeanplot(mod.fit)
autocorr.plot(mod.fit)
geweke.diag(mod.fit)
geweke.plot(mod.fit)
gelman.diag(mod.fit)
heidel.diag(mod.fit)
raftery.diag(mod.fit)
effectiveSize(mod.fit)

### -----caterpillar plot -----------------------------------------
mcmc_intervals(model.mcmc, pars=pars[1:20], prob = 0.8,
               inner_size = 0.9, outer_size = 0.8, prob_outer = 0.95)

#### -----regions meeting criteria ---------------------------------

# extract chains
ex <- MCMCchains(bayes1mcmc, params = 'pi')

# Compute P(pi_i < 0.30) for each column
(probs <- apply(ex, 2, function(x) mean(x < 0.30)))

# Find columns where P(pi_i < 0.30) > 0.9
columns_meeting_criteria <- which(probs > 0.9)

list(
  Columns = columns_meeting_criteria,
  Probabilities = probs[columns_meeting_criteria]
)


#############################################################################
## Model 2 uncentered --------------------------------------------------------
############

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
    beta[j] ~ dnorm(0, 0.001) # age effects
  }
  
  for (k in 1:n_sex) {
    gamma[k] ~ dnorm(0, 0.001) # sex effects
  }
  
  # Binomial likelihood
  for (i in 1:n_munic) {        # Loop over municipalities
    for (j in 1:n_age) {        # Loop over age groups
      for (k in 1:n_sex) {      # Loop over genders
  Yiag[i,j,k] ~ dbin(pi[i,j,k], Niag[i,j,k])
  logit(pi[i,j,k]) <- alpha[n_munic_arr[i, j, k]] + beta[j]*age1[i,j,k] + 
  beta[j]*age2[i,j,k] + beta[j]*age3[i,j,k] + beta[j]*age4[i,j,k] + 
  gamma[k]*male[i,j,k]
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
    beta = rnorm(n_age, 1, 1),
    gamma = rnorm(n_sex, 1, 1),
    mu.int = rnorm(1, 0, 1)
  )
}

params <- c("beta", "gamma", "mu.int", "sigma.int", "pi", "alpha")


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

##-------convergence checks  --------------------------------------------
traceplot(mod2.fit)
rmeanplot(mod2.fit)
autocorr.plot(mod2.fit)
geweke.diag(mod2.fit)
geweke.plot(mod2.fit)
gelman.diag(mod2.fit)
heidel.diag(mod2.fit)
raftery.diag(mod2.fit)
effectiveSize(mod2.fit)


#############################################################################
## Model 2 centered --------------------------------------------------------
############

datajags23 <- list(
  n_munic = length(unique(DDKdata$naam)),
  #n_munic_arr = array(1:n_munic, dim = c(300, 5, 2)),
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

n_munic <- length(unique(DDKdata$naam))
n_age <- length(unique(DDKdata$age))
n_sex <- length(unique(DDKdata$sex))


sink("model2cent.txt")
cat("
model {
# Priors
  for (i in 1:n_munic){
    for (j in 1:n_age) {
      for (k in 1:n_sex) {
        # process model
        Yiag[i,j,k] ~ dbin(pi[i,j,k], Niag[i,j,k])
        logit(pi[i,j,k]) <- b[i,j,k]
        
        # data model
        b[i,j,k] ~ dnorm(mu[i,j,k], tau)
        mu[i,j,k] <- alpha + beta[j]*age1[i,j,k] + 
                          beta[j]*age2[i,j,k] + beta[j]*age3[i,j,k] +
                          beta[j]*age4[i,j,k] + gamma[k]*male[i,j,k]
      }
    }
  }
  
  # priors for fixed effects
  alpha ~ dnorm(0, 0.001) # intercept
  for (j in 1:n_age) {
    beta[j] ~ dnorm(0, 0.001) # Age effects
  }
  
  for (k in 1:n_sex) {
    gamma[k] ~ dnorm(0, 0.001) # Sex effects
  }
  
  tau ~ dgamma(0.001, 0.001) # variance of random effects
}  
", fill = TRUE)
sink()

inits2c <- function() {
  list(
    alpha = rnorm(1, 0, 2),         # Single intercept parameter
    beta = rnorm(n_age, 0, 1),      # Age effects (vector of size n_age)
    gamma = rnorm(n_sex, 0, 1),     # Sex effects (vector of size n_sex)
    tau = rgamma(1, 0.001, 0.001)   # Precision for random intercepts
  )
}

params2c <- c("alpha", "beta", "gamma", "tau")


## Run the model -- r2jags
mod2c.fit <- jags(
  data = datajags23,
  #inits = inits2c,
  parameters.to.save = params2c,
  model.file = "model2cent.txt",
  n.chains = 3,
  n.iter = 10000,
  n.burnin = 5000,
  jags.seed = 123,
  quiet = FALSE
)

# print pD and DIC ---------------------------------------
options(max.print=999999)
sink("model2cresults.txt")
print(mod2c.fit)
sink()

###------ using rjags 
mod2cent.fit <- jags.model(file="model2cent.txt",
                           data = datajags23,
                           n.chains = 3,
                           n.adapt=10000)

update(mod2cent.fit,5000)

model2c.sim <- coda.samples(model = mod2cent.fit,
                            variable.names = params2c,
                            n.iter=10000, 
                            thin=1)

options(max.print=999999)
sink("model2centres.txt")
summary(model2c.sim)
sink()


model2centered.mcmc <- as.mcmc.list(model2c.sim)

##-------convergence checks  --------------------------------------------
traceplot(mod2c.fit)
rmeanplot(mod2c.fit)
autocorr.plot(mod2c.fit)
geweke.diag(mod2c.fit)
geweke.plot(mod2c.fit)
gelman.diag(mod2c.fit)
heidel.diag(mod2c.fit)
raftery.diag(mod2c.fit)
effectiveSize(mod2c.fit) 


########################################################################
### PART 2
########################################################################





