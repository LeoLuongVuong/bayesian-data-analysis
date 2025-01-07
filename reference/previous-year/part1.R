library(rjags)
library(readxl)
library(coda)
library(ggmcmc)
library(ggplot2)


# Part 1
## Model 1
data <- read_excel("D:/bayesian2/new/DDK - 2022.xlsx")
invited <- grep("Invited", colnames(data))
participant <- grep("Participant", colnames(data))
N <- as.matrix(data[, invited])
Y <- as.matrix(data[, participant])
N_array <- array(0, dim = c(num_of_municipalities, AgeGroup, Gender))
Y_array <- array(0, dim = c(num_of_municipalities, AgeGroup, Gender))
for (i in 1:num_of_municipalities) {
  for (a in 1:AgeGroup) {
    for (g in 1:Gender) {
      col_idx <- (a - 1) * Gender + g
      N_array[i, a, g] <- N[i, col_idx]
      Y_array[i, a, g] <- Y[i, col_idx]
    }
  }
}
num_of_municipalities <- nrow(data)
AgeGroup <- 5
Gender <- 2

cat("model {
  for (i in 1:num_of_municipalities) {
    for (a in 1:AgeGroup) {
      for (g in 1:Gender) {
        Y[i, a, g] ~ dbin(pi[i], N[i, a, g])
      }
    }
    pi[i] ~ dbeta(1, 1)
  }
}
", file="model1.txt")

data_new <- list(
  Y = Y_array,
  N = N_array,
  num_of_municipalities= num_of_municipalities,
  AgeGroup = AgeGroup,
  Gender = Gender
)

my.inits <- list(
  list(pi = rep(0.5, num_of_municipalities)),
  list(pi = rep(0.5, num_of_municipalities)),
  list(pi = rep(0.5, num_of_municipalities))
)

parameters <- c("pi")


jags1 <- jags.model(file="model1.txt",
                    data = data_new,
                    inits = my.inits,
                    n.chains = 3)
update(jags1,5000)
model1.sim <- coda.samples(model = jags1,
                           variable.names = parameters,
                           n.iter=10000, 
                           thin=1)

### Model Evaluation
model1.mcmc <- as.mcmc.list(model1.sim)
summary1 <- summary(model1.mcmc)
#traceplot(model1.mcmc[, 1:10])
effectiveSize(model1.mcmc) #(High ESS (e.g., > 1000) suggests good convergence and reliable posterior estimates)
gelman.diag(model1.mcmc) # (target: < 1.1)
geweke.diag(model1.mcmc) #(no convergece)


### Caterpillar Plot
out.ggs1<-ggs(model1.mcmc)
ggs_caterpillar(out.ggs1, model_labels = data$NAAM)

###  P (πi < 0.30|Y ) > 0.9
posterior_pi <- as.matrix(model1.sim)[, grepl("pi", colnames(as.matrix(model1.sim)))]
p <- apply(posterior_pi, 2, function(x) mean(x < 0.30))
result <- p > 0.9
df <- data.frame(Municipality = data$NAAM,
                 Probability = p,
                 Prob_greater_than_0.9 = result)
result <- df[df$Prob_greater_than_0.9 == TRUE, ]
result


## Model 2
cat("
model {
  # Likelihood for Y[iag]
  for (i in 1:num_of_municipalities) {
    for (a in 1:AgeGroup) {
      for (g in 1:Gender) {
        Y[i, a, g] ~ dbin(pi[i, a, g], N[i, a, g])
        logit(pi[i, a, g]) <- alpha + beta * (a-1) + gamma * (g-1) + b[i]
      }
    }
  }

  # Priors for fixed effects
  alpha ~ dnorm(0, 0.01)   # Vague prior for intercept
  beta ~ dnorm(0, 0.01)    # Vague prior for the effect of AgeGroup
  gamma ~ dnorm(0, 0.01)   # Vague prior for the effect of Gender
 

  # Priors for the random effects
  for (i in 1:num_of_municipalities) {
    b[i] ~ dnorm(0, sigma_sq)
  }

  # Prior for the variance of random effects
  sigma_sq ~ dgamma(0.1, 0.1)
}
", file="model2_uncentered.txt")

parameters <- c("pi", "alpha", "beta", "gamma", "b", "sigma_sq")
jags2_uncentered <- jags.model(file="model2_uncentered.txt",
                               data = data_new,
                               n.chains = 3)
update(jags2_uncentered,5000)
model2_uncentered.sim <- coda.samples(model = jags2_uncentered,
                                      variable.names = parameters,
                                      n.iter=10000, 
                                      thin=1)
model2_uncentered.mcmc <- as.mcmc.list(model2_uncentered.sim)
summary2_uncentered <- summary(model2_uncentered.mcmc)


cat("
model {
  # Likelihood for Y[iag]
  for (i in 1:num_of_municipalities) {
    for (a in 1:AgeGroup) {
      for (g in 1:Gender) {
        Y[i, a, g] ~ dbin(pi[i, a, g], N[i, a, g])
        logit(pi[i, a, g]) <- b[i, a, g]
      }
    }
  }

  for (i in 1:num_of_municipalities) {
    for (a in 1:AgeGroup) {
      for (g in 1:Gender) {
        b[i, a, g] ~ dnorm(mu[i, a, g], sigma_sq)
        mu[i, a, g] <- alpha + beta * (a - 1) + gamma * (g - 1)
      }
    }
  }

  # Priors for fixed effects
  alpha ~ dnorm(0, 0.01)   # Vague prior for intercept
  beta ~ dnorm(0, 0.01)    # Vague prior for the effect of AgeGroup
  gamma ~ dnorm(0, 0.01)   # Vague prior for the effect of Gender

  # Prior for the variance of random effects
  sigma_sq ~ dgamma(0.1, 0.1)  # Vague prior for variance of b[iag]
}
", file="model2_centered.txt")

jags2_centered <- jags.model(file="model2_centered.txt",
                             data = data_new,
                             n.chains = 3)
update(jags2_centered,5000)
model2_centered.sim <- coda.samples(model = jags2_centered,
                                    variable.names = parameters,
                                    n.iter=10000, 
                                    thin=1)
model2_centered.mcmc <- as.mcmc.list(model2_centered.sim)
summary2_centered <- summary(model2_centered.mcmc)


### Convergence Checking
traceplot(model2_uncentered.mcmc)
traceplot(model2_centered.mcmc)
effectiveSize(model2_uncentered.mcmc) 
effectiveSize(model2_centered.mcmc) 
gelman.diag(model2_uncentered.mcmc) 
gelman.diag(model2_centered.mcmc) 
geweke.diag(model2_uncentered.mcmc) 
geweke.diag(model2_centered.mcmc) 


### P (πi < 0.30|Y ) > 0.9
# uncentered
alpha_samples_uncentered <- as.matrix(model2_uncentered.sim[, "alpha"])
beta_samples_uncentered <- as.matrix(model2_uncentered.sim[, grep("beta", varnames(model2_uncentered.sim))])
gamma_samples_uncentered <- as.matrix(model2_uncentered.sim[, grep("gamma", varnames(model2_uncentered.sim))])
b_samples_uncentered <- as.matrix(model2_uncentered.sim[, grep("b", varnames(model2_uncentered.sim))])

num_samples_uncentered <- nrow(alpha_samples_uncentered)
pi_samples_uncentered <- array(NA, dim = c(num_samples_uncentered, num_of_municipalities, AgeGroup, Gender))

for (s in 1:num_samples_uncentered) {
  for (i in 1:num_of_municipalities) {
    for (a in 1:AgeGroup) {
      for (g in 1:Gender) {
        logit_pi_uncentered <- alpha_samples_uncentered[s] + beta_samples_uncentered[s]*(a-1) + gamma_samples_uncentered[s]*(g-1) + b_samples[s, i]
        pi_samples_uncentered[s, i, a, g] <- 1 / (1 + exp(-logit_pi_uncentered)) # Logistic function
      }
    }
  }
}

p_below_30_uncentered <- apply(pi_samples_uncentered, c(2, 3, 4), function(x) mean(x < 0.30))
flagged_uncentered <- which(p_below_30_uncentered > 0.9, arr.ind = TRUE)
age_groups <- c("50-54", "55-59", "60-64", "65-69", "70-74")
genders <- c("Male", "Female")
result_table_uncentered <- data.frame(
  Municipality = data$NAAM[flagged_uncentered [, 1]],
  AgeGroup = age_groups[flagged_uncentered [, 2]],
  Gender = genders[flagged_uncentered [, 3]],
  Probability = p_below_30_uncentered[flagged_uncentered ]
)
print(result_table_uncentered)

# centered
alpha_samples_centered <- as.matrix(model2_centered.sim[, "alpha"])
beta_samples_centered <- as.matrix(model2_centered.sim[, grep("beta", varnames(model2_centered.sim))])
gamma_samples_centered <- as.matrix(model2_centered.sim[, grep("gamma", varnames(model2_centered.sim))])
b_samples_centered <- as.matrix(model2_centered.sim[, grep("b", varnames(model2_centered.sim))])

num_samples_centered <- nrow(alpha_samples_centered)
pi_samples_centered <- array(NA, dim = c(num_samples_centered, num_of_municipalities, AgeGroup, Gender))

for (s in 1:num_samples_centered) {
  for (i in 1:num_of_municipalities) {
    for (a in 1:AgeGroup) {
      for (g in 1:Gender) {
        logit_pi_centered <- alpha_samples_centered[s] + beta_samples_centered[s]*(a-1) + gamma_samples_centered[s]*(g-1) + b_samples_centered[s, i]
        pi_samples_centered[s, i, a, g] <- 1 / (1 + exp(-logit_pi_centered)) # Logistic function
      }
    }
  }
}

p_below_30_centered <- apply(pi_samples_centered, c(2, 3, 4), function(x) mean(x < 0.30))
flagged_centered <- which(p_below_30_centered > 0.9, arr.ind = TRUE)
result_table_centered <- data.frame(
  Municipality = data$NAAM[flagged_centered[, 1]],
  AgeGroup = age_groups[flagged_centered[, 2]],
  Gender = genders[flagged_centered[, 3]],
  Probability = p_below_30_centered[flagged_centered]
)
print(result_table_centered)

### Shrinkage
raw_estimates <- Y / N
# uncentered
posterior_pi_uncentered <- model2_uncentered.sim[, grep("^pi", varnames(model2_uncentered.sim))]
posterior_means_uncentered <- apply(as.matrix(posterior_pi_uncentered), 2, mean)
posterior_means_uncentered <- array(posterior_means_uncentered, dim = c(num_of_municipalities, AgeGroup, Gender))


comparison_uncentered <- data.frame(
  Municipality = rep(1:num_of_municipalities, each = AgeGroup * Gender),
  AgeGroup = rep(1:AgeGroup, times = num_of_municipalities * Gender),
  Gender = rep(1:Gender, each = AgeGroup, times = num_of_municipalities),
  RawEstimate = c(raw_estimates),
  PosteriorMean = c(posterior_means_uncentered)
)

ggplot(comparison_uncentered, aes(x = RawEstimate, y = PosteriorMean)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(title = "Raw Estimates vs Posterior Means", x = "Raw Estimate", y = "Posterior Mean")
sd_raw_uncentered <- apply(raw_estimates, 1, sd, na.rm = TRUE)
sd_posterior_uncentered <- apply(posterior_means_uncentered, 1, sd, na.rm = TRUE)

shrinkage_uncentered <- data.frame(
  Municipality = 1:num_of_municipalities,
  SD_Raw = sd_raw_uncentered,
  SD_Posterior = sd_posterior_uncentered,
  Shrinkage = sd_raw_uncentered - sd_posterior_uncentered
)

print(shrinkage_uncentered)

# centered
posterior_pi_centered <- model2_centered.sim[, grep("^pi", varnames(model2_centered.sim))]
posterior_means_centered <- apply(as.matrix(posterior_pi_centered), 2, mean)
posterior_means_centered <- array(posterior_means_centered, dim = c(num_of_municipalities, AgeGroup, Gender))


comparison_centered <- data.frame(
  Municipality = rep(1:num_of_municipalities, each = AgeGroup * Gender),
  AgeGroup = rep(1:AgeGroup, times = num_of_municipalities * Gender),
  Gender = rep(1:Gender, each = AgeGroup, times = num_of_municipalities),
  RawEstimate = c(raw_estimates),
  PosteriorMean = c(posterior_means_centered)
)

ggplot(comparison_centered, aes(x = RawEstimate, y = PosteriorMean)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(title = "Raw Estimates vs Posterior Means", x = "Raw Estimate", y = "Posterior Mean")
sd_raw_centered <- apply(raw_estimates, 1, sd, na.rm = TRUE)
sd_posterior_centered <- apply(posterior_means_centered, 1, sd, na.rm = TRUE)

shrinkage_centered <- data.frame(
  Municipality = 1:num_of_municipalities,
  SD_Raw = sd_raw_centered,
  SD_Posterior = sd_posterior_centered,
  Shrinkage = sd_raw_centered - sd_posterior_centered
)

print(shrinkage_centered)
