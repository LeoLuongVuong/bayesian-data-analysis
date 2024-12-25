# Illustration of the use rjags

library("rjags")
library("readr")

# Change working directory according to your settings

setwd("/Users/lucp1330/Library/CloudStorage/GoogleDrive-christel.faes@uhasselt.be/Mijn Drive/ONDERWIJS/Bayesian Inference II/Course material/Code")

# Read in data
data<-list(
  observed= c(
    5,13,18,5,10,18,29,10,15,22,4,11,10,22,13,14,17,21,25,6,11,21,13,5,19,18,14,17,3,10,
    7,3,12,11,6,16,13,6,9,10,4,9,11,12,23,18,12,7,13,12,12,13,6,14,7,18,13,9,6,8,7,6,16,4,6,12,5,5,
    17,5,7,2,9,7,6,12,13,17,5,5,6,12,10,16,10,16,15,18,6,12,6,8,33,15,14,18,25,14,2,73,13,14,6,20,8,
    12,10,3,11,3,11,13,11,13,10,5,18,10,23,5,9,2,11,9,11,6,11,5,19,15,4,8,9,6,4,4,2,12,12,11,9,7,7,
    8,12,11,23,7,16,46,9,18,12,13,14,14,3,9,15,6,13,13,12,8,11,5,9,8,22,9,2,10,6,10,12,9,11,32,5,11,
    9,11,11,0,9,3,11,11,11,5,4,8,9,30,110),
  expected = c(
    6.17,8.44,7.23,5.62,4.18,29.35,11.79,12.35,7.28,9.40,3.77,3.41,8.70,9.57,8.18,4.35,
    4.91,10.66,16.99,2.94,3.07,5.50,6.47,4.85,9.85,6.95,5.74,5.70,2.22,3.46,4.40,4.05,5.74,6.36,5.13,
    16.99,6.19,5.56,11.69,4.69,6.25,10.84,8.40,13.19,9.25,16.98,8.39,2.86,9.70,12.12,12.94,9.77,
    10.34,5.09,3.29,17.19,5.42,11.39,8.33,4.97,7.14,6.74,17.01,5.80,4.84,12.00,4.50,4.39,16.35,6.02,
    6.42,5.26,4.59,11.86,4.05,5.48,13.13,8.72,2.87,2.13,4.48,5.85,6.67,6.11,5.78,12.31,10.56,10.23,
    2.52,6.22,14.29,5.71,37.93,7.81,9.86,11.61,18.52,12.28,5.41,61.96,8.55,12.07,4.29,19.42,8.25,
    12.90,4.76,5.56,11.11,4.76,10.48,13.13,12.94,14.61,9.26,6.94,16.82,33.49,20.91,5.32,6.77,8.70,
    12.94,16.07,8.87,7.79,14.60,5.10,24.42,17.78,4.04,7.84,9.89,8.45,5.06,4.49,6.25,9.16,12.37,8.40,
    9.57,5.83,9.21,9.64,9.09,12.94,17.42,10.29,7.14,92.50,14.29,15.61,6.00,8.55,15.22,18.42,5.77,
    18.37,13.16,7.69,14.61,15.85,12.77,7.41,14.86,6.94,5.66,9.88,102.16,7.63,5.13,7.58,8.00,12.82,
    18.75,12.33,5.88,64.64,8.62,12.09,11.11,14.10,10.48,7.00,10.23,6.82,15.71,9.65,8.59,8.33,6.06,
    12.31,8.91,50.10,288.00),
  n=195)




# Specification of model

cat("model{

for( i in 1 : n ) {
	
    # Poisson likelihood for observed counts
        observed[i] ~ dpois(lambda[i])   
        lambda[i] <- theta[i]*expected[i]
        	
    # Relative risk
        theta[i] ~ dgamma(alpha,beta)    
        	
    # SMR
        smr[i] <- observed[i]/expected[i]
        	
    # Shrinkage factor
        B[i] <- beta/(beta +  expected[i])     
        	
    # Distribution of future observed counts from region i
         predict[i]  ~ dpois(lambda[i])   

}
	
# Prior distributions for population parameters
    alpha ~ dexp(0.1)
    beta ~  dexp(0.1)
     
    # alpha ~ dunif(0,100)
    # beta ~ dunif(0,100)

# Population mean and population variance
    mtheta <- alpha/beta
    vartheta <- alpha/pow(beta,2)
    sdtheta <- sqrt(vartheta)

# Distribution of future observed counts for a particular expected count, e.g. 100
    theta.new  ~ dgamma(alpha,beta)
    lambda.new <- theta.new*100
    observe.total.new  ~ dpois(lambda.new) 

           
}", file="model.txt")



# Prepare data:

my.data <- list(observed = data$observed,
                expected = data$expected,
                n = data$n)


# Initial parameters

my.inits <- list(
  list(theta = rep(1,data$n), alpha=1, beta=1),
  list(theta = rep(1,data$n), alpha=0.5, beta=0.5),
  list(theta = rep(1,data$n), alpha=1.5, beta=1.5)
)

my.inits <- function()(list(theta = rgamma(data$n,1,1),alpha=rexp(1,1),beta=rexp(1,1)))

# Specify parameters to monitor

parameters <- c("theta","alpha","beta","mtheta","sdtheta")
parameters <- c("theta","alpha","beta","B","predict","theta.new","mtheta","sdtheta")

## Running JAGS:

jags <- jags.model(file="model.txt",
                   data = my.data,
                   inits = my.inits,
                   n.chains = 3)
update(jags,5000) # burn-in period
model.sim <- coda.samples(model = jags,
                          variable.names = parameters,
                          n.iter=10000, 
                          thin=1)


################################################################################
library("coda")

pars<-names(model.sim[[1]][1,])

# Convert model.sim into mcmc.list for processing with CODA package
model.mcmc <- as.mcmc.list(model.sim)

# posterior summary
summary(model.mcmc[,c("alpha","beta","mtheta","sdtheta")])
summary(model.mcmc[,c("theta[1]","theta[2]","theta[3]","theta[4]")])
#summary(model.mcmc[,pars[395:589]])

# traceplot & density plot
plot(model.mcmc[,c("alpha","beta")])
plot(model.mcmc[,c("theta[1]","theta[2]","theta[3]")])


################################################################################
# catterpillar plot
library(ggmcmc)

out.ggs<-ggs(model.mcmc)

ggs_histogram(out.ggs,family ="alpha")

ggs_caterpillar(out.ggs,family ="^theta")




################################################################################
library(MCMCvis)

# posterior summary
MCMCsummary(model.sim,params=c("alpha","beta"),round=2)
MCMCsummary(model.sim,params=c("theta"),round=3, HPD=TRUE)
MCMCsummary(model.sim,params=c("theta\\[([1-9]|[1-2][0-9])\\]"),round=3, ISB=FALSE,exact=FALSE)


# traceplot & density plot
MCMCtrace(model.sim,params=c('alpha','beta'),pdf=FALSE)
MCMCtrace(model.sim,params=c('alpha','beta'),pdf=FALSE,ind=TRUE)

MCMCtrace(model.sim,params=c('theta[1]','theta[2]','theta[3]'),ISB=FALSE,exact=TRUE,pdf=FALSE)
MCMCtrace(model.sim,params=c('theta\\[[1-9]\\]'),ISB=FALSE,exact=FALSE,pdf=FALSE)

# catterpillar plot
MCMCplot(model.sim,params=c("theta\\[([1-9]|[1-3][0-9])\\]"),ISB=FALSE,exact=FALSE)



################################################################################
library(bayesplot)

mcmc_areas(model.mcmc,pars=c("alpha","beta"),prob = 0.8)
mcmc_areas(model.mcmc,pars=c("theta[1]","theta[2]","theta[3]","theta[4]"),prob = 0.8)
mcmc_pairs(model.mcmc,pars=c("theta[1]","theta[2]","theta[3]","theta[4]"))
mcmc_combo(model.mcmc,pars=c("theta[1]","theta[2]","theta[3]","theta[4]"))

mcmc_intervals(model.mcmc,pars=pars[1:10],prob = 0.8)
mcmc_intervals(model.mcmc,pars=pars[395:424],prob = 0.8,point_est="mean")


mcmc_intervals(model.mcmc,pars=pars[c(199:240)],prob = 0.8,point_est="mean")


