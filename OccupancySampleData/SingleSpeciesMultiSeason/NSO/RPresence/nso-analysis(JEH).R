# Single Species, MultiSeason Occupancy analyais
rm(list=ls()); setwd('~/../Downloads/reproblemwithrpresencemodeldo_4/')


# Northern Spotted Owl (Strix occidentalis caurina) in California.
# s=55$ sites visited up to K=$ times per season between 1997 and 2001 ($Y=5$).
# Detection probabilities relatively constant within years, but likely different among years.


#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

#  RPresence package

library(readxl)
library(RPresence)
library(ggplot2)

# get the data read in
# Data for detections should be a data frame with rows corresponding to sites
# and columns to visits.
# The usual 1=detected; 0=not detected; NA=not visited is used.

input.history <- read.csv(file.path("..","NSO.csv"), header=FALSE, skip=2, na.strings="-")
input.history$V1 <- NULL # drop the site number

# do some basic checks on your data
# e.g. check number of sites; number of visits etc
nrow(input.history)
ncol(input.history)
range(input.history, na.rm=TRUE)
sum(is.na(input.history))


# Five years with 8 visits. Don't need same number of visits/season
Nvisits.per.season  <- rep(8,5) # five years with 8 visits. Don't need same number of visits/season

# Create the *.pao file
nso.pao <- RPresence::createPao(input.history,
                                nsurveyseason=Nvisits.per.season,
                                title='NSO SSMS')
nso.pao

#-----
# Fit some models.
# We start with the first formulation in terms of psi, gamma, epsilon, p (do.1_)
# Note that formula DO NOT HAVE AN = SIGN

  mod.psiDot.gDot.eDot.pDot <- RPresence::occMod(model=list(psi~1, gamma~1, epsilon~1, p~1),
                                        type="do.1", data=nso.pao)
summary(mod.psiDot.gDot.eDot.pDot)

names(mod.psiDot.gDot.eDot.pDot)
names(mod.psiDot.gDot.eDot.pDot$real)

# Estimate of initial occupance
mod.psiDot.gDot.eDot.pDot$real$psi[1,]

# Estimate of  local extinction probability for each unit
mod.psiDot.gDot.eDot.pDot$real$epsilon[ seq(1, by=nrow(input.history), length.out=length(Nvisits.per.season)-1),]

# Estimate of  local colonization probability for each unit
mod.psiDot.gDot.eDot.pDot$real$gamma[ seq(1, by=nrow(input.history), length.out=length(Nvisits.per.season)-1),]

# Estimate of probability of detection at each time point for each unit
mod.psiDot.gDot.eDot.pDot$real$p[ grepl('unit1_', row.names(mod.psiDot.gDot.eDot.pDot$real$p)),]

# Derived parameters - estimated occupancy for each unit in years 2....
names(mod.psiDot.gDot.eDot.pDot$derived)
mod.psiDot.gDot.eDot.pDot$derived$psi[ grepl('unit1_', row.names(mod.psiDot.gDot.eDot.pDot$derived$psi)),]

# Not possible to obtain lambda (change in occupancy) or lambda' (change in odds of occupancy) directly
# at this point. But we can estimate it. Not simple to get the standard errors, except by bootstrapping everything.
# Stack together all of the psi for a site
site_psi <- rbind(
  mod.psiDot.gDot.eDot.pDot$real$psi   [grepl('unit1_', row.names(mod.psiDot.gDot.eDot.pDot$real   $psi)),],
  mod.psiDot.gDot.eDot.pDot$derived$psi[grepl('unit1_', row.names(mod.psiDot.gDot.eDot.pDot$derived$psi)),])
site_psi

logit <- function(x) log(x/(1-x))
expit <- function(x) 1/(1+exp(-x))

lambda <- exp(diff(log(site_psi$est),1))
lambda
prod(lambda) # overall growth rate over entire set of seasons

lambda.prime <- exp(diff(logit(site_psi$est),1))
lambda.prime
prod(lambda.prime) # overall growth rate on logit scale over entire set of season


#-----
# psi(.) gamma(.) epsilon(.) p(season)
# We start with the first formulation in terms of psi, gamma, epsilon, p (do.1_)
# Note that formula DO NOT HAVE AN = SIGN
# Notice the use of the reserved keyword SEASON
mod.psiDot.gDot.eDot.pYear <- RPresence::occMod(model=list(psi~1, gamma~1, epsilon~1, p~SEASON),
                                        type="do.1", data=nso.pao)
summary(mod.psiDot.gDot.eDot.pYear)

# Estimate of initial occupance
mod.psiDot.gDot.eDot.pYear$real$psi[1,]

# Estimate of  local extinction probability for each unit
mod.psiDot.gDot.eDot.pYear$real$epsilon[ seq(1, by=nrow(input.history), length.out=length(Nvisits.per.season)-1),]

# Estimate of  local colonization probability for each unit
mod.psiDot.gDot.eDot.pYear$real$gamma[ seq(1, by=nrow(input.history), length.out=length(Nvisits.per.season)-1),]

# Estimate of probability of detection at each time point for each unit
mod.psiDot.gDot.eDot.pYear$real$p[ grepl('unit1_', row.names(mod.psiDot.gDot.eDot.pYear$real$p)),]

# Derived parameters - estimated occupancy for each unit in years 2....
names(mod.psiDot.gDot.eDot.pYear$derived)
mod.psiDot.gDot.eDot.pYear$derived$psi[ grepl('unit1_', row.names(mod.psiDot.gDot.eDot.pYear$derived$psi)),]

# Not possible to obtain lambda (change in occupancy) or lambda' (change in odds of occupancy) directly
# at this point. But we can estimate it. Not simple to get the standard errors, except by bootstrapping everything.
# Stack together all of the psi for a site
site_psi <- rbind(
  mod.psiDot.gDot.eDot.pYear$real$psi   [grepl('unit1_', row.names(mod.psiDot.gDot.eDot.pYear$real   $psi)),],
  mod.psiDot.gDot.eDot.pYear$derived$psi[grepl('unit1_', row.names(mod.psiDot.gDot.eDot.pYear$derived$psi)),])
site_psi

logit <- function(x) log(x/(1-x))
expit <- function(x) 1/(1+exp(-x))

lambda <- exp(diff(log(site_psi$est),1))
lambda
prod(lambda) # overall growth rate over entire set of seasons

lambda.prime <- exp(diff(logit(site_psi$est),1))
lambda.prime
prod(lambda.prime) # overall growth rate on logit scale over entire set of season

#-----
# psi(.) gamma(season) epsilon(season) p(season)
# We start with the first formulation in terms of psi, gamma, epsilon, p (do.1_)
# Note that formula DO NOT HAVE AN = SIGN
# Notice the use of the reserved keyword SEASON
mod.psiDot.gYear.eYear.pYear <- RPresence::occMod(model=list(psi~1, gamma~SEASON, epsilon~SEASON, p~SEASON),
                                        type="do.1", data=nso.pao)
summary(mod.psiDot.gYear.eYear.pYear)

# Estimate of initial occupance
mod.psiDot.gYear.eYear.pYear$real$psi[1,]

# Estimate of  local extinction probability for each unit
mod.psiDot.gYear.eYear.pYear$real$epsilon[ seq(1, by=nrow(input.history), length.out=length(Nvisits.per.season)-1),]

# Estimate of  local colonization probability for each unit
mod.psiDot.gYear.eYear.pYear$real$gamma[ seq(1, by=nrow(input.history), length.out=length(Nvisits.per.season)-1),]

# Estimate of probability of detection at each time point for each unit
mod.psiDot.gYear.eYear.pYear$real$p[ grepl('unit1_', row.names(mod.psiDot.gYear.eYear.pYear$real$p)),]

# Derived parameters - estimated occupancy for each unit in years 2....
names(mod.psiDot.gYear.eYear.pYear$derived)
mod.psiDot.gYear.eYear.pYear$derived$psi[ grepl('unit1_', row.names(mod.psiDot.gYear.eYear.pYear$derived$psi)),]

# Not possible to obtain lambda (change in occupancy) or lambda' (change in odds of occupancy) directly
# at this point. But we can estimate it. Not simple to get the standard errors, except by bootstrapping everything.
# Stack together all of the psi for a site
site_psi <- rbind(
  mod.psiDot.gYear.eYear.pYear$real$psi   [grepl('unit1_', row.names(mod.psiDot.gYear.eYear.pYear$real   $psi)),],
  mod.psiDot.gYear.eYear.pYear$derived$psi[grepl('unit1_', row.names(mod.psiDot.gYear.eYear.pYear$derived$psi)),])
site_psi

logit <- function(x) log(x/(1-x))
expit <- function(x) 1/(1+exp(-x))

lambda <- exp(diff(log(site_psi$est),1))
lambda
prod(lambda) # overall growth rate over entire set of seasons

lambda.prime <- exp(diff(logit(site_psi$est),1))
lambda.prime
prod(lambda.prime) # overall growth rate on logit scale over entire set of season

#-----
# Model with no change in occupancy.
# We start with the first formulation in terms of psi, gamma, epsilon, p (do.1_)
# Note that formula DO NOT HAVE AN = SIGN
# Notice the use of the reserved keyword SEASON

# this doesn't seem to work in this version of Presence

# You need to look at the design matrices to know the parameter names (see dmat object)
mod.psiDot.gYear.eYear.pYear$dmat$epsilon

fixed.parm <- data.frame(param=c("gamma1","epsilon1"),     value=0)
fixed.parm
gamcov=data.frame(SEASN=as.factor(c(rep(1,2*55),rep(2:3,each=55))))

mod.psiDot.g0.e0.pYear <- RPresence::occMod(
        model=list(psi~1, gamma~SEASN, epsilon~1, p~SEASON),
        type="do.1", cov.list=list(gamma.cov=gamcov),
        fixed=fixed.parm,
        data=nso.pao,outfile='modname');

summary(mod.psiDot.g0.e0.pYear)

mod.psiDot.g0.e0.pYear$dmat
# Estimate of initial occupance
mod.psiDot.g0.e0.pYear$real$psi[1,]

# Estimate of  local extinction probability for each unit
mod.psiDot.g0.e0.pYear$real$epsilon[ seq(1, by=nrow(input.history), length.out=length(Nvisits.per.season)-1),]

# Estimate of  local colonization probability for each unit
mod.psiDot.g0.e0.pYear$real$gamma[ seq(1, by=nrow(input.history), length.out=length(Nvisits.per.season)-1),]

# Estimate of probability of detection at each time point for each unit
mod.psiDot.g0.e0.pYear$real$p[ grepl('unit1_', row.names(mod.psiDot.g0.e0.pYear$real$p)),]

# Derived parameters - estimated occupancy for each unit in years 2....
names(mod.psiDot.g0.e0.pYear$derived)
mod.psiDot.g0.e0.pYear$derived$psi[ grepl('unit1_', row.names(mod.psiDot.g0.e0.pYear$derived$psi)),]

#-----
# Random occupancy are fit using type="do.4" in the call.
# Parameters are psi,  p with gamma=1-epsilon enforced internally
# Note that formula DO NOT HAVE AN = SIGN

mod.psiYear.pYear.RO <- RPresence::occMod(
        model=list(psi~SEASON, p~SEASON),
        type="do.4",
        data=nso.pao,outfile='modname')

summary(mod.psiYear.pYear.RO)

# Estimate of occupancy
mod.psiYear.pYear.RO$real$psi[ grepl('unit1_', row.names(mod.psiYear.pYear.RO$real$psi)),]

# Estimate of probability of detection at each time point for each unit
mod.psiYear.pYear.RO$real$p[ grepl('unit1_', row.names(mod.psiYear.pYear.RO$real$p)),]

# For some reason these are not given in the current version or RPresence.
# Derived parameters - estimated occupancy for each unit in years 2....
names(mod.psiYear.pYear.RO$derived)
mod.psiYear.pYear.RO$derived$psi[ grepl('unit1_', row.names(mod.psiYear.pYear.RO$derived$psi)),]

# Derived parameter Estimate of  local extinction probability for each unit
mod.psiYear.pYear.RO$derived$epsilon[ seq(1, by=nrow(input.history), length.out=length(Nvisits.per.season)-1),]

# Estimate of  local colonization probability for each unit
mod.psiYear.pYear.RO$real$gamma[ seq(1, by=nrow(input.history), length.out=length(Nvisits.per.season)-1),]

#------
# Model averaging
models<-list(mod.psiDot.gDot.eDot.pDot,
             mod.psiDot.gDot.eDot.pYear,
             mod.psiDot.gYear.eYear.pYear,
             mod.psiYear.pYear.RO)
results<-RPresence::createAicTable(models)
summary(results)

nunits=models[[1]]$data$nunits
nseasons=models[[1]]$data$nseasons

# unable to model average other year psi values, or lambda values (groan) unless you write your own
cat('psi:\n')
print(RPresence::modAvg(results, param="psi")[1:5,])

models<-list(mod.psiDot.gDot.eDot.pDot,
             mod.psiDot.gDot.eDot.pYear,
             mod.psiDot.gYear.eYear.pYear)
results<-RPresence::createAicTable(models)
i=seq(1, nunits*(nseasons-1), nunits)
print(RPresence::modAvg(results, param="epsilon")[i,])

# ================== compute model-averaged estimates of derived psi's... ==================
est=var=wgts=NULL

for (i in 1:length(models)) {
  est=cbind(est,models[[i]]$derived$psi[,"est"])   #  estimate of psi from each model
  var=cbind(var,models[[i]]$derived$psi[,"se"]^2)  #  variance of psi from each model
  j=which(models[[i]]$modname==results$table$Model) #  get model weight of each model
  wgts=c(wgts,results$table$wgt[j])
}
modavgest=rowSums(t(wgts*t(est)))  #  model-avg estimate = sum (wgt(i)*estimate(i))
var2=(est-modavgest)^2                              #  model-avg variance = sum(wgt(i)*sqrt(var(i) + (est(i)-mod.avg.estimate)^2))
modavgse= sapply(1:nrow(est),function(i) sum(wgts*sqrt(var[i,] + var2[i,])))       #  model-avg se = sqrt(mod.avg.var)
modavgpsis=cbind(modavgest,modavgse)
rownames(modavgpsis)=paste0('psi',rep(2:nseasons,each=nunits),'_unit',rep(1:nunits,nseasons-2))
i=grep('unit1$',rownames(modavgpsis))
print(modavgpsis[i,])
#  ==========================================================================================

#  ========  compute model-averaged estimates of lambda(1)  (psi(2)/psi(1))  ================
#               (must do separately since psi1 is a real estimted parameter and psi2 is a derived parameter)
modavgout=NULL
#                 get model weight for each model...
wgts=sapply(1:length(models),function(i) results$table$wgt[which(models[[i]]$modname==results$table$Model)])

for (j in 1:nunits) {   #  for each site, compute lam(1) = psi(2)/psi(1)
  est=var=NULL
  for (i in 1:length(models)) {              #  for each model...
    psi1=models[[i]]$real$psi[j,"est"]          #  get psi(1) for site j
    psi2=models[[i]]$derived$psi[j,"est"]       #  get psi(2) for site j
    lam=psi2/psi1                               #  compute lambda = psi2/psi1
    #                              get varianc-covariance matrix of psi1,psi2
    vcm=diag(c(models[[i]]$real$psi[j,"se"]^2,models[[i]]$derived$psi[j,"se"]^2))

    grad=c(1/psi1,-psi2/psi1^2)      #  compute derivative if lambda with respect to psi1,psi2
    vari=t(grad) %*% vcm %*% grad    #  delta method to compute variance of lambda
    est=c(est,lam)                   #  save estimate of lambda for each model
    var=c(var,vari)                  #  save estimate of variance for each model
  }
  maest=sum(wgts * est)                       #  model-avg estimate = sum (wgt(i)*estimate(i))
  var2=(est-maest)^2                          #  model-avg variance = sum(wgt(i)*sqrt(var(i) + (est(i)-mod.avg.estimate)^2))
  mase= sum(wgts*sqrt(var + var2))       #  model-avg se = sqrt(mod.avg.var)
  modavgout=rbind(modavgout,c(maest,mase))
}
#  ==========================================================================================

#  ========  compute model-averaged estimates of lambda(1)  (psi(i+1)/psi(i) for i>1 ================

for (k in 2:(nseasons-1)) {
  for (j in 1:nunits) {   #  for each site, compute lam(1) = psi(2)/psi(1)
    est=var=NULL
    for (i in 1:length(models)) {              #  for each model...
      psi1=models[[i]]$derived$psi[(k-2)*nunits+j,"est"]       #  get psi(i) for site j
      psi2=models[[i]]$derived$psi[(k-1)*nunits+j,"est"]       #  get psi(i+1) for site j
      lam=psi2/psi1                               #  compute lambda = psi2/psi1
      #                              get varianc-covariance matrix of psi1,psi2
      vcm=diag(c(models[[i]]$derived$psi[(k-2)*nunits+j,"se"]^2,models[[i]]$derived$psi[(k-1)*nunits+j,"se"]^2))

      grad=c(1/psi1,-psi2/psi1^2)      #  compute derivative if lambda with respect to psi1,psi2
      vari=t(grad) %*% vcm %*% grad    #  delta method to compute variance of lambda
      est=c(est,lam)                   #  save estimate of lambda for each model
      var=c(var,vari)                  #  save estimate of variance for each model
    }
    maest=sum(wgts * est)                       #  model-avg estimate = sum (wgt(i)*estimate(i))
    var2=(est-maest)^2                          #  model-avg variance = sum(var(i)*wgt(i) + wgt(i)*(est(i)-mod.avg.estimate)^2)
    mase= sqrt(sum(var*wgts + var2*wgts))       #  model-avg se = sqrt(mod.avg.var)
    modavgout=rbind(modavgout,c(maest,mase))
  }
}
rownames(modavgout)=paste0('lambda',rep(1:(nseasons-1),each=nunits),'_unit',rep(1:nunits,nseasons-1))
colnames(modavgout)=c('estimate','std.err')
i=grep('unit1$',rownames(modavgout))         #  only print for unit1... change this if you want others...
print(modavgout[i,])

