
# Illustration of the use rjags

library("rjags")
library("coda")
library("readr")
library("runjags")

# Change working directory according to your settings

setwd("/Users/lucp1330/Library/CloudStorage/GoogleDrive-christel.faes@uhasselt.be/Mijn Drive/ONDERWIJS/Bayesian Inference II/Course material/Code")

# Read in data
source("ToenailData1.R")

# Visualisation of data
par(mfrow=c(1,2))
plot(data$time[data$treat==0],data$response[data$treat==0],main="intraconazol",xlab="time",ylab="Response",ylim=c(0,24))
u<-unique(data$id[data$treat==0])
for (i in 1:length(u)){lines(data$time[data$id==u[i]],data$response[data$id==u[i]],col="grey")}
for (i in 1:3){lines(data$time[data$id==u[i]],data$response[data$id==u[i]],col="blue")}
m0<-NULL
for (i in 1:length(unique(data$time))){
  m0[i]<-mean(data$response[data$treat==0 & data$time==(unique(data$time))[i]])
}
lines(unique(data$time),m0,col="red",lwd=3)

plot(data$time[data$treat==1],data$response[data$treat==1],main="lamisil",xlab="time",ylab="Response",ylim=c(0,24))
u<-unique(data$id[data$treat==1])
for (i in 1:length(u)){lines(data$time[data$id==u[i]],data$response[data$id==u[i]],col="grey")}
for (i in 1:3){lines(data$time[data$id==u[i]],data$response[data$id==u[i]],col="blue")}
m1<-NULL
for (i in 1:length(unique(data$time))){
  m1[i]<-mean(data$response[data$treat==1 & data$time==(unique(data$time))[i]])
}
lines(unique(data$time),m1,col="red",lwd=3)


# Specification of model

cat("model{
  

  # level 1
  
    for( iobs in 1 : Nobs ) {	
     	mean[iobs] <- beta0 + beta1*time[iobs] + beta2*treat[iobs]*time[iobs]+ 
	                          b[id[iobs],1] + b[id[iobs],2]*time[iobs]
	    response[iobs] ~ dnorm(mean[iobs], tau.eps)
	
    # prediction new observation from current set of observations
    
      newresp[iobs] ~ dnorm(mean[iobs], tau.eps)
	}
	
	# level 2
	
    for( isubj in 1 : Nsubj ) {	
	    	b[isubj,1:2] ~ dmnorm(meanb[1:2],taub[1:2,1:2])
	  }
			
  	# prediction of new observation but not from current set 
	 
	  bnew[1:2] ~ dmnorm(meanb[1:2],taub[1:2,1:2])
	  for( iobs in 1 : Nnew ) {	
	      meannew[iobs] <- beta0 + beta1*timenew[iobs] + beta2*treatnew[iobs]*timenew[iobs]+ 
	                          bnew[1] + bnew[2]*timenew[iobs]
	       nnewresp[iobs] ~ dnorm(meannew[iobs], tau.eps)
	      }
	
 # prior distributions
	 
    taub[1:2,1:2] ~ dwish(Rb[1:2,1:2],2)
    sigma2b[1:2,1:2]  <- inverse(taub[1:2,1:2])
		for (i in 1:2) {
	    	sigmab[i] <- sqrt(sigma2b[i,i])
	  }
		corrb <- sigma2b[1,2]/(sigmab[1]*sigmab[2])
		
	  tau.eps ~ dgamma(0.001,0.001)
	  sigma.eps <- 1 / sqrt(tau.eps)
	  sigma.eps2 <- pow(sigma.eps,2)
		
	  beta0  ~ dnorm(0.0,1.0E-6)
	  beta1  ~ dnorm(0.0,1.0E-6)
	  beta2  ~ dnorm(0.0,1.0E-6)
	  
  }", file="model.txt")


#file.show("model.txt")


# Prepare data:

my.data <- list(Nobs=data$Nobs, Nsubj=data$Nsubj,
                time=data$time, treat=data$treat,
                response=data$response,id=data$id,
                meanb=data$meanb, Rb=data$Rb,
                Nnew=data$Nnew, timenew=data$timenew, treatnew=data$treatnew)


# Initial parameters
my.inits <- list(
  list(beta0 = 0, beta1 = 0, beta2 = 0, tau.eps=1),
  list(beta0 = 0, beta1 = 0, beta2 = 0, tau.eps=10),
  list(beta0 = 0, beta1 = 0, beta2 = 0, tau.eps=5)
)

# Specify parameters to monitor

parameters <- c("beta0","beta1","beta2","sigma.eps","corrb","sigmab","newresp","nnewresp")

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
#traceplot(model.mcmc)

# BGR diagnostic (target: < 1.1)

#gelman.plot(model.mcmc,ask=FALSE)
gelman.diag(model.mcmc)

# Geweke diagnostic

geweke.diag(model.mcmc)
#geweke.plot(osteo.mcmc,ask=FALSE)



# Plot of posterior predictive distribution
par(mfrow=c(2,2))
densplot(model.mcmc[,"newresp[1]"])
abline(v=data$response[1],col="red")  
densplot(model.mcmc[,"newresp[2]"])
abline(v=data$response[2],col="red")  
densplot(model.mcmc[,"newresp[3]"])
abline(v=data$response[3],col="red")  
densplot(model.mcmc[,"newresp[4]"])
abline(v=data$response[4],col="red")  


par(mfrow=c(2,2))
densplot(model.mcmc[,"nnewresp[1]"])
abline(v=m0[1],col="red")  
densplot(model.mcmc[,"nnewresp[2]"])
abline(v=m0[2],col="red")  
densplot(model.mcmc[,"nnewresp[3]"])
abline(v=m0[3],col="red")  
densplot(model.mcmc[,"nnewresp[4]"])
abline(v=m0[4],col="red")  

