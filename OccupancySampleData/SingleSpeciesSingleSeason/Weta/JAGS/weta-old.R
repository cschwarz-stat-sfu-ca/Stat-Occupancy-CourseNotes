
#################################################################################
## R-code for the Bayesian Mahoenui Giant Weta example in MacKenzie et al. (2018).
## Example illustrates the effect of different prior distributions on the analysis
## and results. Data is read in from an exisiting PRESENCE data file, weta.pao.
## User is expected to understand R, OpenBUGS/JAGS code.
#################################################################################

## Set working directory. Will need to be edited for your computer to where
## data file is stored.
setwd("Weta")

# Load required R libraries
library(RPresence) # not available on CRAN. https://www.mbr-pwrc.usgs.gov/software/presence.shtml
library(jagsUI)
library(R2OpenBUGS)

# Read in data file
weta<-readPao("weta.pao")

# Extract detection data and covariates from 'weta' object.
h<-as.matrix(weta$det.data) # detection data

browse<-as.numeric(weta$unitcov$Browsed) # browse covariate, =1 if plot was browsed
obs2<-weta$survcov$Obs2 # Obs2 covariate, =1 if survey conducted by observer 2
obs3<-weta$survcov$Obs3 # Obs3 covariate, =1 if survey conducted by observer 3

# Set any missing covariate values to 0. Specific value doesn't matter as no
# detection data associated with that survey.
obs2[is.na(obs2)]<-0 
obs3[is.na(obs3)]<-0

# Reformat observer covariates to a 'wide' format, the same as the detection data.
obs2<-matrix(obs2,ncol=5)
obs3<-matrix(obs3,ncol=5)

#####################################################################
## Define JAGS model with vague priors for regression coefficients,
## variance = 100, or precision = 1/variance = 0.01.
#####################################################################
model<-function(){
  
  ## model structure
  for (ii in 1:s){
    logit(psi[ii])<-a1+a2*browse[ii]
    z[ii]~dbern(psi[ii])
    for(jj in 1:K){
      logit(p[ii,jj])<-b1+b2*obs2[ii,jj]+b3*obs3[ii,jj]
      p.temp[ii,jj]<-p[ii,jj]*z[ii]
      h[ii,jj]~dbern(p.temp[ii,jj])
    }
  }
  
  a1~dnorm(0,0.01)
  a2~dnorm(0,0.01)
  
  b1~dnorm(0,0.01)
  b2~dnorm(0,0.01)
  b3~dnorm(0,0.01)
  
  logit(psi1[1])<-a1
  logit(psi1[2])<-a1+a2
  
  logit(p1[1])<-b1
  logit(p1[2])<-b1+b2
  logit(p1[3])<-b1+b3
}

write.model(model,con="model_vague.txt") # write JAGS model code to file
mcmc.data<-list(h=h,s=72,K=5, browse=browse,obs2=obs2,obs3=obs3) # define input variables for model
mcmc.params<-c("psi1","p1","a1","a2","b1","b2","b3") # define which node values to store
mcmc.inits<-function() {list(z=rep(1,72))}

Sys.time()
vague <- jags(data=mcmc.data, inits=mcmc.inits, parameters.to.save=mcmc.params,
              n.iter=11000, model.file="model_vague.txt",n.chains=3,parallel=TRUE,verbose=TRUE,n.burnin = 1000)
Sys.time()

#####################################################################
## Define JAGS model with semi-informative priors for regression coefficients,
## variance = 10, or precision = 1/variance = 0.1.
#####################################################################
model<-function(){
  
  ## model structure
  for (ii in 1:s){
    logit(psi[ii])<-a1+a2*browse[ii]
    z[ii]~dbern(psi[ii])
    for(jj in 1:K){
      logit(p[ii,jj])<-b1+b2*obs2[ii,jj]+b3*obs3[ii,jj]
      p.temp[ii,jj]<-p[ii,jj]*z[ii]
      h[ii,jj]~dbern(p.temp[ii,jj])
    }
  }
  
  a1~dnorm(0,0.1)
  a2~dnorm(0,0.1)

  b1~dnorm(0,0.1)
  b2~dnorm(0,0.1)
  b3~dnorm(0,0.1)
  
  logit(psi1[1])<-a1
  logit(psi1[2])<-a1+a2
  
  logit(p1[1])<-b1
  logit(p1[2])<-b1+b2
  logit(p1[3])<-b1+b3
  
}

write.model(model,con="model_vague.txt") # write JAGS model code to file
mcmc.data<-list(h=h,s=72,K=5, browse=browse,obs2=obs2,obs3=obs3) # define input variables for model
mcmc.params<-c("psi1","p1","a1","a2","b1","b2","b3") # define which node values to store
mcmc.inits<-function() {list(z=rep(1,72))}


write.model(model,con="model_semi.txt") # write JAGS model code to file
mcmc.data<-list(h=h,s=72,K=5, browse=browse,obs2=obs2,obs3=obs3) # define input variables for model
mcmc.params<-c("psi1","p1","a1","a2","b1","b2","b3") # define which node values to store
mcmc.inits<-function() {list(z=rep(1,72))}

Sys.time()
semi <- jags(data=mcmc.data, inits=mcmc.inits, parameters.to.save=mcmc.params,
                n.iter=11000, model.file="model_semi.txt",n.chains=3,parallel=TRUE,verbose=TRUE,n.burnin = 1000)
Sys.time()

#####################################################################
## Define JAGS model with informative priors for regression coefficients,
## variance = 1, or precision = 1/variance = 1.
#####################################################################
model<-function(){
  
  ## model structure
  for (ii in 1:s){
    logit(psi[ii])<-a1+a2*browse[ii]
    z[ii]~dbern(psi[ii])
    for(jj in 1:K){
      logit(p[ii,jj])<-b1+b2*obs2[ii,jj]+b3*obs3[ii,jj]
      p.temp[ii,jj]<-p[ii,jj]*z[ii]
      h[ii,jj]~dbern(p.temp[ii,jj])
    }
  }
  
  a1~dnorm(0,1)
  a2~dnorm(0,1)
  
  b1~dnorm(0,1)
  b2~dnorm(0,1)
  b3~dnorm(0,1)
  
  logit(psi1[1])<-a1
  logit(psi1[2])<-a1+a2
  
  logit(p1[1])<-b1
  logit(p1[2])<-b1+b2
  logit(p1[3])<-b1+b3
  
}

write.model(model,con="model_inf.txt") # write JAGS model code to file
mcmc.data<-list(h=h,s=72,K=5, browse=browse,obs2=obs2,obs3=obs3) # define input variables for model
mcmc.params<-c("psi1","p1","a1","a2","b1","b2","b3") # define which node values to store
mcmc.inits<-function() {list(z=rep(1,72))}

Sys.time()
informed <- jags(data=mcmc.data, inits=mcmc.inits, parameters.to.save=mcmc.params,
              n.iter=11000, model.file="model_inf.txt",n.chains=3,parallel=TRUE,verbose=TRUE,n.burnin = 1000)
Sys.time()

##########################################################################################
## Numerical summaries of posterior distribtuions
##########################################################################################
vague
semi
informed

##########################################################################################
## Example plots
##########################################################################################
jpeg("weta_trace%1d.jpg",width=1200,height=600,res=144,quality=100)
traceplot(vague$samples[,"a2"],ylab=expression(alpha[2]))
traceplot(vague$samples[,"psi1[2]"],ylab=expression(psi[Browse]),ylim=c(0,1))
dev.off()

#####################################################################
## Define JAGS model with informative priors for regression coefficients,
## variance = 1, or precision = 1/variance = 1.
## Standardise covariates
#####################################################################
browse<-browse-mean(browse)
obs2<-obs2-mean(obs2,na.rm=TRUE)
obs3<-obs3-mean(obs3,na.rm=TRUE)

model<-function(){
  
  ## model structure
  for (ii in 1:s){
    logit(psi[ii])<-a1+a2*browse[ii]
    z[ii]~dbern(psi[ii])
    for(jj in 1:K){
      logit(p[ii,jj])<-b1+b2*obs2[ii,jj]+b3*obs3[ii,jj]
      p.temp[ii,jj]<-p[ii,jj]*z[ii]
      h[ii,jj]~dbern(p.temp[ii,jj])
    }
  }
  
  a1~dt(0,0.16,3)
  a2~dt(0,0.16,3)
  
  b1~dt(0,0.16,3)
  b2~dt(0,0.16,3)
  b3~dt(0,0.16,3)

  logit(psi1[1])<-a1+a2*-0.4861111
  logit(psi1[2])<-a1+a2*0.5138889 
  
  logit(p1[1])<-b1+b2*-0.3091603+b3*-0.3435115
  logit(p1[2])<-b1+b2*0.6908397+b3*-0.3435115
  logit(p1[3])<-b1+b2*-0.3091603+b3*0.6564885
  
}

write.model(model,con="model_gelman.txt") # write JAGS model code to file
mcmc.data<-list(h=h,s=72,K=5, browse=browse,obs2=obs2,obs3=obs3) # define input variables for model
mcmc.params<-c("psi1","p1","a1","a2","b1","b2","b3") # define which node values to store
mcmc.inits<-function() {list(z=rep(1,72))}

Sys.time()
gelman <- jags(data=mcmc.data, inits=mcmc.inits, parameters.to.save=mcmc.params,
               n.iter=11000, model.file="model_gelman.txt",n.chains=3,parallel=TRUE,verbose=TRUE,n.burnin = 1000)
Sys.time()

gelman
