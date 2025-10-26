# single species multi season 

# Northern Spotted Owl (Strix occidentalis caurina) in California.
# s=55 sites visited up to K=4 times per season between 1997 and 2001 Y=5).
# Detection probabilities relatively constant within years, but likely different among years.

# Using the unmarked package

# 2018-08-17 Code contributed by Carl James Schwarz (cschwarz.stat.sfu.ca@gmail.com)

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

# unmarked package

library(readxl)
library(unmarked)
library(ggplot2)

# get the data read in
# Data for detections should be a data frame with rows corresponding to sites
# and columns to visits.
# The usual 1=detected; 0=not detected; NA=not visited is used.
# Each year must have the same number of visit. Pad with missing if this is not true.

input.history <- read.csv(file.path("..","NSO.csv"), header=FALSE, skip=2, na.strings="-")
input.history$V1 <- NULL # drop the site number

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.history)
ncol(input.history)
range(input.history, na.rm=TRUE)
sum(is.na(input.history))

# Five years with 8 visits. input.history must have same number of visits/year
Nsites = nrow(input.history)
Nyears = 5
Nvisits= ncol(input.history)/Nyears


# Create the site covariates
site.covar <- data.frame(SiteNum=1:Nsites)
head(site.covar)

# Create the Yearly site (site x time) covariates. 
# yearlySiteCovs is a data frame with Nsites x Nyears rows which are in site-major, year-minor order.
# I.e. list covariates for Site 1 for years 1..., then Site 2 for years .... etc
yearsite.covar <- data.frame(SiteNum=rep(1:Nsites,  each=Nyears),
                             Year   =as.factor(rep(1:Nyears,Nsites)))
head(yearsite.covar)

# Create the Observation covarates (site x time x visit)
# Observation covaraites is a data frame with Nsite x Ntime x Nvisit rows
# which are ordered by site-year-observation, so that a column of obsCovs corresponds to as.vector(t(y))
test.history <- matrix(1:40, nrow=4, ncol=10, byrow=TRUE)
test.history
as.vector(t(test.history))

obs.covar <- data.frame( SiteNum = rep(1:Nsites, each=Nyears * Nvisits),
                         Year   =  as.factor(rep( rep(1:Nyears, each=Nvisits), Nsites)),
                         Visit  =  as.factor(rep(1:Nvisits, Nyears*Nsites)))
obs.covar[1:30,]

# Create the UMF file
nso.UMF <- unmarked::unmarkedMultFrame(input.history,
                                siteCovs=site.covar,
                                yearlySiteCovs=yearsite.covar,
                                obsCovs = obs.covar,
                                numPrimary=Nyears)
nso.UMF
summary(nso.UMF)
plot(nso.UMF)


#-----
# Fit some models.
# We start with the most basic model for psi, gamma, epsilon, p 
# Note that formula DO  HAVE AN = SIGN now
mod.psiDot.gDot.eDot.pDot <- unmarked::colext(
     psiformula= ~1, 
     gammaformula = ~ 1, 
     epsilonformula = ~ 1,
     pformula = ~ 1,
     data=nso.UMF,
     se=TRUE)
summary(mod.psiDot.gDot.eDot.pDot)

slotNames(mod.psiDot.gDot.eDot.pDot)


# We get some estimates using predict in the usual way

names(mod.psiDot.gDot.eDot.pDot)
newdata <- data.frame(intercept=1)
predict(mod.psiDot.gDot.eDot.pDot, type='psi', newdata=newdata )

predict(mod.psiDot.gDot.eDot.pDot, type='col', newdata=newdata )

predict(mod.psiDot.gDot.eDot.pDot, type='ext', newdata=newdata )

predict(mod.psiDot.gDot.eDot.pDot, type='det', newdata=newdata )

# estimates of psi moving ahead.
# There are two types.
#   - population level occupancy, i.e. for the sampled population from which the sites were selected
#   - actual site level occupancy, i.e. for the sites actually sampled.

# What is the occupancy for the population of ALL potential sites
# first index = not occupied/occupied, second index=year, last index=site
mod.psiDot.gDot.eDot.pDot@projected[,,1]
mod.psiDot.gDot.eDot.pDot@projected[2,,1]  
projected(mod.psiDot.gDot.eDot.pDot)  # this gives the mean population occupancy over all sites

# What is the actual occupancy for these particular sites
mod.psiDot.gDot.eDot.pDot@smoothed[2,,1:2]
smoothed(mod.psiDot.gDot.eDot.pDot)[2,] # this gives the mean actual occupancy over all sites


# Obtain standard errors for these forecasts using parametric bootstrapping
psi.projected <- function(model, site) {
  logit <- function(x) log(x/(1-x))
  expit <- function(x) 1/(1+exp(-x))
  
  psi.projected <- model@projected[2,,site]
  
  # compute population growth
  lambda <- exp(diff(log(psi.projected),1))
  lambda.overall <- prod(lambda) # overall growth rate over entire set of seasons
  
  lambda.prime <- exp(diff(logit(psi.projected),1))
  lambda.prime.overall <- prod(lambda.prime) # overall growth rate on logit scale over entire set of season
  c(psi.projected=psi.projected, 
    lambda=       lambda, 
    lambda.overall=lambda.overall, 
    lambda.prime  =lambda.prime, 
    lambda.prime.overall=lambda.prime.overall )
}
psi.projected(mod.psiDot.gDot.eDot.pDot, site=1)

# use bootstrapping to find standard errors.
# You likely want to do 1000 bootstrap samples
mod.psiDot.gDot.eDot.pDot.boot <- parboot(mod.psiDot.gDot.eDot.pDot, statistic=psi.projected, nsim=100, site=1)
slotNames(mod.psiDot.gDot.eDot.pDot.boot)
mod.psiDot.gDot.eDot.pDot.boot@t0
cbind(est=mod.psiDot.gDot.eDot.pDot.boot@t0,
      se=apply(mod.psiDot.gDot.eDot.pDot.boot@t.star, 2, sd),
      t(apply(mod.psiDot.gDot.eDot.pDot.boot@t.star, 2, quantile, probs=c(0.025, 0.975))))


#-----
# Fit some models.
# Model for ps(.)i, gamma(.), epsilon(.), p(year) 
# Note that formula DO  HAVE AN = SIGN now
mod.psiDot.gDot.eDot.pYear <- unmarked::colext(
  psiformula= ~1, 
  gammaformula = ~ 1, 
  epsilonformula = ~ 1,
  pformula = ~ Year,
  data=nso.UMF,
  se=TRUE)
summary(mod.psiDot.gDot.eDot.pYear)

slotNames(mod.psiDot.gDot.eDot.pYear)


# We get some estimates using predict in the usual way

names(mod.psiDot.gDot.eDot.pYear)
newdata <- data.frame(Year=as.factor(1:5))
predict(mod.psiDot.gDot.eDot.pYear, type='psi', newdata=newdata )

predict(mod.psiDot.gDot.eDot.pYear, type='col', newdata=newdata )

predict(mod.psiDot.gDot.eDot.pYear, type='ext', newdata=newdata )


predict(mod.psiDot.gDot.eDot.pYear, type='det', newdata=newdata )

# estimates of psi moving ahead.
# There are two types.
#   - population level occupancy, i.e. for the sampled population from which the sites were selected
#   - actual site level occupancy, i.e. for the sites actually sampled.

# What is the occupancy for the population of ALL potential sites
# first index = not occupied/occupied, second index=year, last index=site
mod.psiDot.gDot.eDot.pYear@projected[,,1]
mod.psiDot.gDot.eDot.pYear@projected[2,,1]  
projected(mod.psiDot.gDot.eDot.pYear)  # this gives the mean population occupancy over all sites

# What is the actual occupancy for these particular sites
mod.psiDot.gDot.eDot.pYear@smoothed[2,,1:2]
smoothed(mod.psiDot.gDot.eDot.pYear)[2,] # this gives the mean actual occupancy over all sites


# Obtain standard errors for these forecasts using parametric bootstrapping
psi.projected <- function(model, site) {
  logit <- function(x) log(x/(1-x))
  expit <- function(x) 1/(1+exp(-x))
  
  psi.projected <- model@projected[2,,site]
  
  # compute population growth
  lambda <- exp(diff(log(psi.projected),1))
  lambda.overall <- prod(lambda) # overall growth rate over entire set of seasons
  
  lambda.prime <- exp(diff(logit(psi.projected),1))
  lambda.prime.overall <- prod(lambda.prime) # overall growth rate on logit scale over entire set of season
  c(psi.projected=psi.projected, 
    lambda=       lambda, 
    lambda.overall=lambda.overall, 
    lambda.prime  =lambda.prime, 
    lambda.prime.overall=lambda.prime.overall )
}
psi.projected(mod.psiDot.gDot.eDot.pYear, site=1)

# use bootstrapping to find standard errors.
# You likely want to do 1000 bootstrap samples
mod.psiDot.gDot.eDot.pYear.boot <- parboot(mod.psiDot.gDot.eDot.pYear, statistic=psi.projected, nsim=100, site=1)
slotNames(mod.psiDot.gDot.eDot.pYear.boot)
mod.psiDot.gDot.eDot.pYear.boot@t0
cbind(est=mod.psiDot.gDot.eDot.pYear.boot@t0,
      se=apply(mod.psiDot.gDot.eDot.pYear.boot@t.star, 2, sd),
      t(apply(mod.psiDot.gDot.eDot.pYear.boot@t.star, 2, quantile, probs=c(0.025, 0.975))))


#-----
# Fit some models.
# Model for ps(.)i, gamma(Year), epsilon(Year), p(year) 
# Note that formula DO  HAVE AN = SIGN now
mod.psiDot.gYear.eYear.pYear <- unmarked::colext(
  psiformula= ~1, 
  gammaformula = ~ Year, 
  epsilonformula = ~ Year,
  pformula = ~ Year,
  data=nso.UMF,
  se=TRUE)
summary(mod.psiDot.gYear.eYear.pYear)

slotNames(mod.psiDot.gYear.eYear.pYear)


# We get some estimates using predict in the usual way

names(mod.psiDot.gYear.eYear.pYear)
newdata <- data.frame(Year=as.factor(1:5))
predict(mod.psiDot.gYear.eYear.pYear, type='psi', newdata=newdata )
predict(mod.psiDot.gYear.eYear.pYear, type='det', newdata=newdata )

newdata <- data.frame(Year=as.factor(2:5))
predict(mod.psiDot.gYear.eYear.pYear, type='col', newdata=newdata )

predict(mod.psiDot.gYear.eYear.pYear, type='ext', newdata=newdata )


# estimates of psi moving ahead.
# There are two types.
#   - population level occupancy, i.e. for the sampled population from which the sites were selected
#   - actual site level occupancy, i.e. for the sites actually sampled.

# What is the occupancy for the population of ALL potential sites
# first index = not occupied/occupied, second index=year, last index=site
mod.psiDot.gYear.eYear.pYear@projected[,,1]
mod.psiDot.gYear.eYear.pYear@projected[2,,1]  
projected(mod.psiDot.gYear.eYear.pYear)  # this gives the mean population occupancy over all sites

# What is the actual occupancy for these particular sites
mod.psiDot.gYear.eYear.pYear@smoothed[2,,1:2]
smoothed(mod.psiDot.gYear.eYear.pYear)[2,] # this gives the mean actual occupancy over all sites


# Obtain standard errors for these forecasts using parametric bootstrapping
psi.projected <- function(model, site) {
  logit <- function(x) log(x/(1-x))
  expit <- function(x) 1/(1+exp(-x))
  
  psi.projected <- model@projected[2,,site]
  
  # compute population growth
  lambda <- exp(diff(log(psi.projected),1))
  lambda.overall <- prod(lambda) # overall growth rate over entire set of seasons
  
  lambda.prime <- exp(diff(logit(psi.projected),1))
  lambda.prime.overall <- prod(lambda.prime) # overall growth rate on logit scale over entire set of season
  c(psi.projected=psi.projected, 
    lambda=       lambda, 
    lambda.overall=lambda.overall, 
    lambda.prime  =lambda.prime, 
    lambda.prime.overall=lambda.prime.overall )
}
psi.projected(mod.psiDot.gYear.eYear.pYear, site=1)

# use bootstrapping to find standard errors.
# You likely want to do 1000 bootstrap samples
mod.psiDot.gYear.eYear.pYear.boot <- parboot(mod.psiDot.gYear.eYear.pYear, statistic=psi.projected, nsim=100, site=1)
slotNames(mod.psiDot.gYear.eYear.pYear.boot)
mod.psiDot.gYear.eYear.pYear.boot@t0
cbind(est=mod.psiDot.gYear.eYear.pYear.boot@t0,
      se=apply(mod.psiDot.gYear.eYear.pYear.boot@t.star, 2, sd),
      t(apply(mod.psiDot.gYear.eYear.pYear.boot@t.star, 2, quantile, probs=c(0.025, 0.975))))



#------
# Model averaging
models.u <-unmarked::fitList(
  mod.psiDot.gDot.eDot.pDot,
  mod.psiDot.gDot.eDot.pYear,
  mod.psiDot.gYear.eYear.pYear)
aic.table.u <- unmarked::modSel(models.u)
aic.table.u

# Get model averaged estimates of occupancy
predict(models.u, type="psi")[1,]

# Get model averaged estimates of detection. 
# Notice these are one big list of nsites x nvisits long
# with the detection probabilities for each site
cbind(obs.covar, predict(models.u, type="det"))[1:10,]

# Not possible to model average projected psi's directly.
# You need to do this on your own (groan)


