#library(rjags)
library(R2jags)
library(coda)
library(readxl)
library(pacman)
library(MCMCvis)
library(mcmcplots)

pacman::p_load(tidyverse, coda,readxl,janitor)

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

# write.csv(bda_data_long, 'data/bda_data_long', row.names = FALSE)


###################################################################

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

# Niag = array(c(bda_data_long$invited), dim = c(300, 5, 2))
# Niag
# n_munic <- length(unique(bda_data_long$naam))
# n_age <- length(unique(bda_data_long$age))
# n_sex <- length(unique(bda_data_long$sex))

sink()

sink("model1.txt")
cat("
# ...
model {
  for (i in 1:n_munic) {        # Loop over municipalities
    pi[i] ~ dbeta(1, 1)         # Prior for participation rate per municipality
    
    for (j in 1:n_age) {        # Loop over age groups
      for (k in 1:n_sex) {      # Loop over genders
        Yiag[i,j,k] ~ dbin(pi[i], Niag[i,j,k])  # Likelihood
      }
    }
  }
}
", fill = TRUE)

sink()


# inits for chains
# inits1 <- list(pi = rep(0.1, 300))
# inits2 <- list(pi = rep(0.5, 300))
# 
# inits <- function() {
#   list(pi = rep(0.1, 300))
# }

# Initial parameters
my.inits <- list(
  list(pi = rep(0.1, 300)),
  list(pi = rep(0.5, 300)),
  list(pi = rep(0.9, 300))
)

# params to monitor
params <- c("pi")

# Run the model
mod.fit <- jags(
  data = datajags,
  inits = my.inits, 
  parameters.to.save = c("pi"),
  model.file = "model1.txt",
  n.chains = 3,
  n.iter = 2000,
  n.burnin = 1000,
  jags.seed = 123,
  quiet = FALSE
)

traceplot(mod.fit)

bayes.mod.fit.mcmc <- as.mcmc(mod.fit)
summary(bayes.mod.fit.mcmc)

autocorr.plot(bayes.mod.fit.mcmc)

gelman.plot(bayes.mod.fit.mcmc)
geweke.diag(bayes.mod.fit.mcmc)

heidel.diag(bayes.mod.fit.mcmc)

caterplot(bayes.mod.fit.mcmc, parms = c("pi"))
# Leo: this looks good but I will try to make it a bit prettier 

# caterplot(bayes.mod.fit.mcmc, parms = c("pi[1]","pi[2]","pi[3]"))

# Leo: General comment: I would definitely try to finetune the code for the MCMC diagnostics (i.e. lines 110 to 122)

### MCMCvis
MCMCsummary(bayes.mod.fit.mcmc, round = 3)

MCMCplot(bayes.mod.fit.mcmc, 
         params = 'pi',
         rank = TRUE,
         guide_lines = TRUE)

MCMCtrace(bayes.mod.fit.mcmc, 
          params = c('pi[1]', 'pi[2]', 'pi[3]'), 
          ISB = FALSE, 
          exact = TRUE,
          pdf = FALSE)

### Extract posterior samples for pi

# extract chains
ex <- MCMCchains(bayes.mod.fit.mcmc, params = 'pi')

# Compute P(pi_i < 0.30) for each column
(probs <- apply(ex, 2, function(x) mean(x < 0.30)))
# Leo: indeed here you got 0 but it's because of the error you made earlier

# Find columns where P(pi_i < 0.30) > 0.9
columns_meeting_criteria <- which(probs > 0.9)

# Display results
list(
  Columns = columns_meeting_criteria,
  Probabilities = probs[columns_meeting_criteria]
)

