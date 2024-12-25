
# Illustration of the use rjags

library("rjags")
library("coda")
library("readr")
library("runjags")

# Change working directory according to your settings

setwd("/Users/lucp1330/Library/CloudStorage/GoogleDrive-christel.faes@uhasselt.be/Mijn Drive/ONDERWIJS/Bayesian Inference II/Course material/Code")

# Read in data
source("ToenailData2.R")


# Specification of model

cat("model{
	
	# level 1
	
	  for( iobs in 1 : Nobs ) {	

       logit(p[iobs]) <- beta[1] + beta[2]*time[iobs] + beta[3]*treat[iobs]*time[iobs]+ 
	                             b0[id[iobs]]  +  b1[id[iobs]] *time[iobs]

       response[iobs] ~ dbern(p[iobs])
	      }
	
	# level 2
	
	 for( isubj in 1 : Nsubj ) {	
			b0[isubj] ~ dnorm(0,taub0)
			b1[isubj] ~ dnorm(0,taub1)
	 }
			
			
	# priors
	
		sigmab0 ~ dunif(0,100)
		sigma2b0  <- pow(sigmab0,2)
		taub0 <- pow(sigma2b0,-1)
		sigmab1 ~ dunif(0,10)
		sigma2b1  <- pow(sigmab1,2)
		taub1 <- pow(sigma2b1,-1)
		
		for (i in 1:3){
		beta[i]  ~ dnorm(0.0,1.0E-2)
		}
	  
  }", file="model.txt")


#file.show("model.txt")


# Prepare data:

my.data <- list(Nobs=data$Nobs, Nsubj=data$Nsubj,
                time=data$time, treat=data$treat,
                response=data$response2,id=data$id)


# Initial parameters
my.inits <- list(
  list(beta=c(0,0,0) , sigmab0=1, sigmab1=1),
  list(beta=c(0,0,0) , sigmab0=10, sigmab1=1),
  list(beta=c(0,0,0) , sigmab0=1, sigmab1=10)
)

# Specify parameters to monitor

parameters <- c("beta","sigmab0","sigmab1")

## Running JAGS:

jags <- jags.model(file="model.txt",
                   data = my.data,
                   inits = my.inits,
                   n.chains = 3)
update(jags,1000) # burn-in period
model.sim <- coda.samples(model = jags,
                          variable.names = parameters,
                          n.iter=30000, 
                          thin=3)



# Convert osteo.sim into mcmc.list for processing with CODA package

model.mcmc <- as.mcmc.list(model.sim)

# Produce general summary of obtained MCMC sampling

summary(model.mcmc)

# Trace plots from Gibbs sampler

par(mfrow=c(1,2))
traceplot(model.mcmc[,"beta[1]"])
traceplot(model.mcmc[,"beta[2]"])
traceplot(model.mcmc[,"beta[3]"])

# BGR diagnostic (target: < 1.1)

#gelman.plot(model.mcmc,ask=FALSE)
gelman.diag(model.mcmc)

# Geweke diagnostic

geweke.diag(model.mcmc)
#geweke.plot(osteo.mcmc,ask=FALSE)



