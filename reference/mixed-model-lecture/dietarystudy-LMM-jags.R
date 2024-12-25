
# Illustration of the use rjags

library("rjags")
library("coda")
library("readr")
library("runjags")

# Change working directory according to your settings

# setwd("/Users/lucp1330/Library/CloudStorage/GoogleDrive-christel.faes@uhasselt.be/Mijn Drive/ONDERWIJS/Bayesian Inference II/Course material/Code")

# Read in data
dietary.study<-read.table(file = "reference/mixed-model-lecture/cholesterol.txt",header=T)

# Specification of model

cat("model{
  # C = number of centers
	# nc[i] = number of subjects in center i
  # the mean of center i is given in theta[i] of length 8
	# tau is 1/sigma^2
	
	# likelihood - level 1
      for (j in 1:N){
	      chol[j] ~ dnorm(theta[cluster[j]],tau)
	    }
	 
	# the distribution of the thetas - level 2
      for (i in 1:C){
	      theta[i] ~ dnorm(mu,tau.theta) 
	      B[i] <- tau.theta/(tau.theta+ nc[i]*tau)
      }
	    
	# prediction for cluster 1
	    for (j in 1:nc[1]){
	    predict[j] ~ dnorm(theta[cluster[j]],tau)
	    }
	    mpredict<-mean(predict[])
		
  # prior for sigma2
		  tau ~ dgamma(0.001,0.001)
		  sigma <- 1 / sqrt(tau)
		  sigma.2 <- 1/tau
		
  # prior for mu
	  	mu ~ dnorm(0.0,1.0E-6)	 
		
	# prior for sigma.theta
		  sigma.theta ~ dunif(0,100)
		  sigma.theta2 <- pow(sigma.theta,2)
		  tau.theta<- pow(sigma.theta,-2)
		
	# intra-class correlation
		  r <- sigma.theta2/(sigma.theta2+sigma.2)
		
	# prediction new observations of a new subsidiary
		
		  theta.new ~ dnorm(mu,tau.theta)
		  for (j in 1:nnew){  
	      	predict.new[j] ~ dnorm(theta.new,tau)
			 }
		  mpredict.new <- mean(predict.new[])
  }", file="model.txt")


#file.show("model.txt")


# Prepare data:

my.data <- list(chol=dietary.study$chol,
                cluster=dietary.study$cluster,
                N=nrow(dietary.study),
                C=length(unique(dietary.study$cluster)),
                nnew=82,
                nc=c(82,51,71,71,62,69,74,83))


# Initial parameters
C<-length(unique(dietary.study$cluster))
my.inits <- list(
  list(theta = rep(0,C),mu = 0, tau=1, sigma.theta=1),
  list(theta = rep(0,C),mu = 0, tau=10, sigma.theta=1),
  list(theta = rep(0,C),mu = 0, tau=1, sigma.theta=10)
)

# Specify parameters to monitor

parameters <- c("theta","mu","sigma","sigma.theta","mpredict.new","r","B","mpredict")

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
traceplot(model.mcmc)

# BGR diagnostic (target: < 1.1)

#gelman.plot(model.mcmc,ask=FALSE)
gelman.diag(model.mcmc)

# Geweke diagnostic

geweke.diag(model.mcmc)
#geweke.plot(osteo.mcmc,ask=FALSE)

# Plot of posterior distribution
par(mfrow=c(1,1))
densplot(model.mcmc[,"mu"])


# catterpillar plot
library(ggmcmc)

out.ggs<-ggs(model.mcmc)
ggs_caterpillar(out.ggs,family ="^theta",sort=FALSE)
ggs_caterpillar(out.ggs,family ="^theta",sort=TRUE)


# Plot of posterior predictive distribution
par(mfrow=c(1,2))
densplot(model.mcmc[,"mpredict.new"])
abline(v=mean(combine.mcmc(model.mcmc)[,"mu"]),col="red")  

densplot(model.mcmc[,"mpredict"])
abline(v=mean(combine.mcmc(model.mcmc)[,"theta[1]"]),col="red")  

