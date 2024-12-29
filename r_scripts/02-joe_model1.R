#library(rjags)
library(R2jags)
library(coda)
library(readxl)
library(pacman)
library(MCMCvis)
library(mcmcplots)
library(ggmcmc)

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

bda_data_long <- read.csv("data/bda_data_long.csv")

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


sink("model1.txt")
cat("
# ...
model {
  for (i in 1:n_munic) {        # Loop over municipalities
    pi[i] ~ dbeta(0.5, 0.5)         # Prior for participation rate per municipality
    
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
  n.iter = 10000,
  n.burnin = 5000,
  jags.seed = 123,
  quiet = FALSE
)

par("mar")
par(mar=c(1,1,1,1))


traceplot(mod.fit)

bayes1mcmc <- as.mcmc(mod.fit)
summary(bayes1mcmc)

jags <- jags.model(file="model1.txt",
                   data = datajags,
                   inits = my.inits,
                   n.chains = 3)
summary(jags)

update(jags, 5000) # burn-in period

model.sim <- coda.samples(model = jags,
                          variable.names = params,
                          n.iter=10000, 
                          thin=1)
sink("model1res2.txt")
summary(model.sim)
sink()
#######################################################################
library(coda)


pars <- names(model.sim[[1]][1,])
# Convert model.sim into mcmc.list for processing with CODA package
model.mcmc <- as.mcmc.list(model.sim)

# posterior summary
summary(model.mcmc)

# traceplot & density plot
plot(model.mcmc)
plot(model.mcmc, cex.axis = 0.4)

# autocorrelation and running mean plots
png("pictures/autocorr.png", width = 30, 
    height = 20, units = "cm", res = 300)
autocorr.plot(model.sim)
dev.off()
rmeanplot(model.sim)


# geweke diagnostics
geweke.diag(model.sim)
par("mar")
par(mar=c(1,1,1,1))
png("pictures/geweke.png", width = 30, 
    height = 20, units = "cm", res = 300)
geweke.plot(model.sim)
dev.off()

# BGR diagnostic
gelman.diag(model.sim)
png("pictures/GRB.png", width = 30, 
    height = 20, units = "cm", res = 300)
gelman.plot(model.sim)
dev.off()

# HW diagnostic
sink("model1res.txt")
heidel.diag(model.sim)
sink()

# Raftery–Lewis (RL) diagnostic
raftery.diag(model.sim)


# trace and density plots
plot(bayes1mcmc)
ggs_traceplot(bayes1mcmcggs, family = "^pi")


#######################################################################
## GGMCMC

# catterpillar plot
library(ggmcmc)

# caterpillar plot 
bayes1mcmcggs <- ggs(bayes1mcmc)

ggs_traceplot(bayes1mcmcggs, family = "^pi")


ggs_histogram(bayes1mcmcggs,family ="^pi")

png("pictures/caterpillar.png", width = 18, 
    height = 9, units = "cm", res = 300)
ggs_caterpillar(bayes1mcmcggs, family = "^pi")
dev.off()

################################################################################
library(MCMCvis)

# traceplot & density plot
png("pictures/trace_.png", width = 18, 
    height = 9, units = "cm", res = 300)
MCMCtrace(model.sim,params='pi',pdf=FALSE)
dev.off()

# subset of traceplots
png("pictures/trace6-9.png", width = 18, 
    height = 9, units = "cm", res = 300)
MCMCtrace(model.sim, params=c('pi\\[[6-9]\\]'),ISB=FALSE,exact=FALSE,pdf=FALSE)
dev.off()

MCMCsummary(model.sim)
#########################################################################
########## MCMCPLOTS

library(mcmcplots)
par("mar")
par(mar=c(1,1,1,1))

png("pictures/caterpillar.png", width = 50, 
    height = 40, units = "cm", res = 400)
caterplot(model.sim)
dev.off()

rmeanplot(bayes1mcmc)

#######################################################################
## BAYESPLOT
library(bayesplot)

# caterpillar plot
png("pictures/cater1-20.png", width = 18, 
    height = 9, units = "cm", res = 300)
mcmc_intervals(model.mcmc, pars=pars[1:20], prob = 0.8, inner_size = 0.9, outer_size = 0.8, prob_outer = 0.95) +
  labs(
    subtitle = "Caterpillar plot: pi[1:20]"
  )
dev.off()
png("pictures/cater21-40.png", width = 18, 
    height = 9, units = "cm", res = 300)
mcmc_intervals(model.mcmc, pars=pars[21:40], prob = 0.8, inner_size = 0.9, outer_size = 0.8, prob_outer = 0.95) +
  labs(
    subtitle = "Caterpillar plot: pi[21:40]"
  )
dev.off()
png("pictures/cater41-60.png", width = 18, 
    height = 9, units = "cm", res = 300)
mcmc_intervals(model.mcmc, pars=pars[41:60], prob = 0.8, inner_size = 0.9, outer_size = 0.8, prob_outer = 0.95) +
  labs(
    subtitle = "Caterpillar plot: pi[41:60]"
  )
dev.off()
png("pictures/cater61-80.png", width = 18, 
    height = 9, units = "cm", res = 300)
mcmc_intervals(model.mcmc, pars=pars[61:80], prob = 0.8, inner_size = 0.9, outer_size = 0.8, prob_outer = 0.95) +
  labs(
    subtitle = "Caterpillar plot: pi[61:80]"
  )
dev.off()
png("pictures/cater81-100.png", width = 18, 
    height = 9, units = "cm", res = 300)
mcmc_intervals(model.mcmc, pars=pars[81:100], prob = 0.8, inner_size = 0.9, outer_size = 0.8, prob_outer = 0.95) +
  labs(
    subtitle = "Caterpillar plot: pi[81:100]"
  )
dev.off()
png("pictures/cater101-120.png", width = 18, 
    height = 9, units = "cm", res = 300)
mcmc_intervals(model.mcmc, pars=pars[101:120], prob = 0.8, inner_size = 0.9, outer_size = 0.8, prob_outer = 0.95) +
  labs(
    subtitle = "Caterpillar plot: pi[101:120]"
  )
dev.off()
png("pictures/cater121-140.png", width = 18, 
    height = 9, units = "cm", res = 300)
mcmc_intervals(model.mcmc, pars=pars[121:140], prob = 0.8, inner_size = 0.9, outer_size = 0.8, prob_outer = 0.95) +
  labs(
    subtitle = "Caterpillar plot: pi[121:140]"
  )
dev.off()
png("pictures/cater141-160.png", width = 18, 
    height = 9, units = "cm", res = 300)
mcmc_intervals(model.mcmc, pars=pars[141:160], prob = 0.8, inner_size = 0.9, outer_size = 0.8, prob_outer = 0.95) +
  labs(
    subtitle = "Caterpillar plot: pi[141:160]"
  )
dev.off()
png("pictures/cater161-180.png", width = 18, 
    height = 9, units = "cm", res = 300)
mcmc_intervals(model.mcmc, pars=pars[161:180], prob = 0.8, inner_size = 0.9, outer_size = 0.8, prob_outer = 0.95) +
  labs(
    subtitle = "Caterpillar plot: pi[161:180]"
  )
dev.off()
png("pictures/cater181-200.png", width = 18, 
    height = 9, units = "cm", res = 300)
mcmc_intervals(model.mcmc, pars=pars[181:200], prob = 0.8, inner_size = 0.9, outer_size = 0.8, prob_outer = 0.95) +
  labs(
    subtitle = "Caterpillar plot: pi[181:200]"
  )
dev.off()
png("pictures/cater201-220.png", width = 18, 
    height = 9, units = "cm", res = 300)
mcmc_intervals(model.mcmc, pars=pars[201:220], prob = 0.8, inner_size = 0.9, outer_size = 0.8, prob_outer = 0.95) +
  labs(
    subtitle = "Caterpillar plot: pi[201:220]"
  )
dev.off()
png("pictures/cater221-240.png", width = 18, 
    height = 9, units = "cm", res = 300)
mcmc_intervals(model.mcmc, pars=pars[221:240], prob = 0.8, inner_size = 0.9, outer_size = 0.8, prob_outer = 0.95) +
  labs(
    subtitle = "Caterpillar plot: pi[221:240]"
  )
dev.off()
png("pictures/cater241-260.png", width = 18, 
    height = 9, units = "cm", res = 300)
mcmc_intervals(model.mcmc, pars=pars[241:260], prob = 0.8, inner_size = 0.9, outer_size = 0.8, prob_outer = 0.95) +
  labs(
    subtitle = "Caterpillar plot: pi[241:260]"
  )
dev.off()
png("pictures/cater261-280.png", width = 18, 
    height = 9, units = "cm", res = 300)
mcmc_intervals(model.mcmc, pars=pars[261:280], prob = 0.8, inner_size = 0.9, outer_size = 0.8, prob_outer = 0.95) +
  labs(
    subtitle = "Caterpillar plot: pi[261:280]"
  )
dev.off()
png("pictures/cater281-300.png", width = 18, 
    height = 9, units = "cm", res = 300)
mcmc_intervals(model.mcmc, pars=pars[281:300], prob = 0.8, inner_size = 0.9, outer_size = 0.8, prob_outer = 0.95) +
  labs(
    subtitle = "Caterpillar plot: pi[281:300]"
  )
dev.off()

######################################################################

### Extract posterior samples for pi

# extract chains
ex <- MCMCchains(bayes1mcmc, params = 'pi')

# Compute P(pi_i < 0.30) for each column
(probs <- apply(ex, 2, function(x) mean(x < 0.30)))s

# Find columns where P(pi_i < 0.30) > 0.9
columns_meeting_criteria <- which(probs > 0.9)

# Display results
list(
  Columns = columns_meeting_criteria,
  Probabilities = probs[columns_meeting_criteria]
)
#############################################################################
library(dplyr)
library(readxl)

# summary statistics
model1res <- read_excel("model1res.xlsx")

# municipality with max and min rate
(maxpi <- model1res %>%
  filter(mean == min(mean)) %>%
  select(pi, mean))

# Count how many pi have mean >= 0.50
(count <- model1res %>%
  filter(mean >= 0.50) %>%
  summarise(count = n()))

# Count how many pi have mean < 0.25
(count <- model1res %>%
    filter(mean < 0.25) %>%
    #summarise(count = n()) %>%
    select(pi, mean))
