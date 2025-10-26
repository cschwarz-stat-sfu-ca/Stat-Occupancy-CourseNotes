# Single Species Single Season Occupancy models 
# Weta Example
# Using unmarked


# 2018-08-15 Code contributed by Carl James Schwarz (cschwarz.stat.sfu.cs@gmail.com)

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------


library(readxl)
library(unmarked)
library(ggplot2)

# get the data read in
# Data for detections should be a data frame with rows corresponding to sites
# and columns to visits.
# The usual 1=detected; 0=not detected; NA=not visited is used.

input.history <- readxl::read_excel(file.path("..","weta.xls"),
                                 sheet="detection_histories",
                                 na="-",
                                 col_names=FALSE)  # notice no column names in row 1 of data file. 

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.history)
ncol(input.history)
range(input.history, na.rm=TRUE)
sum(is.na(input.history))

head(input.history)

# Get the site level covariates
site_covar <- readxl::read_excel(file.path("..","weta.xls"),
                                 sheet="site_covar",
                                 na="-",
                                 col_names=TRUE)  # notice col_names in row 1 of table. 


# Create an alternate site level covariate that is a categorical variable rather 
# than indicator variables
site_covar$BrowCat <- paste(c("","B")[1+unlist(site_covar[,1])], c("","N")[1+unlist(site_covar[,2])], sep="")
xtabs(~BrowCat, data=site_covar,exclude=NULL, na.action=na.pass)
colSums(site_covar[,1:2])

head(site_covar)

# Visit Covariates
# unmarked does not have any reserved keywords like RPresence so 
# we need to construct a covariate data frame with the covariates
# again stacked. The ordering is different in unmarked.
# Here you want for site 1, the visit specific covariates, then
# for site 2, etc.
# This needs to be added to the unmarked dataframe.

# Get the individual covariates. 
obs1 <- readxl::read_excel(file.path("..","weta.xls"),
                                 sheet="Obs1",
                                 na="-",
                                 col_names=FALSE) 
obs2 <- readxl::read_excel(file.path("..","weta.xls"),
                                 sheet="Obs2",
                                 na="-",
                                 col_names=FALSE) 
obs3 <- readxl::read_excel(file.path("..","weta.xls"),
                                 sheet="Obs3",
                                 na="-",
                                 col_names=FALSE) 

obs <- obs1*1 + obs2*2 + obs3*3
head(obs)


# We also want to set up covariate values where we
# set the detection probability equal in first 2 occasions and last 3 occasions.
# We need to define a survey covariate that has two levels 
survey.covar.u <- data.frame( Site=rep(1:nrow(input.history), each=ncol(input.history)),
                         Visit=rep(1:ncol(input.history), nrow(input.history)),
                         Time =rep(c("T1","T2","T3","T4","T5"), nrow(input.history)),
                         obs1 =as.vector(t(as.matrix(obs1))),
                         obs2 =as.vector(t(as.matrix(obs2))),
                         obs3 =as.vector(t(as.matrix(obs3))),
                         obs  =as.character(as.vector(t(as.matrix(obs)))), stringsAsFactors=FALSE)
                         
survey.covar.u[1:10,]

# check that missing values in history and observer covariates align
select <- is.na(as.vector(t(as.matrix(input.history))))
survey.covar.u[select,]
sum(is.na(survey.covar.u[!select,]))


weta.UMF <- unmarked::unmarkedFrameOccu(
                     y = input.history,
                     siteCovs=site_covar,
                     obsCovs=survey.covar.u
                     ) 
summary(weta.UMF)
plot(weta.UMF)



#-----
# Fit some models.
# Note that formula DO NOT HAVE AN = SIGN
# The two formula are for detecton and occupancy in that order
mod.pdot.u <- unmarked::occu(~1 ~1 , data=weta.UMF, se=TRUE)

# This creates a complicated S4 object with SLOTS that can be accessed using various
# functions. This is more complcated than with RPresence
slotNames(mod.pdot.u)
mod.pdot.u@estimates



# These functions return estimates on the LOGIT scale - not very useful for most people
coef(mod.pdot.u, type="state")
SE(mod.pdot.u, type="state")
confint(mod.pdot.u, type="state")

# It is possible to use backTransform() but only if NO covariates
backTransform(mod.pdot.u, type='state')

# better to use predict() to do the estimation.
# Because there are no covariates here, the data.frame is empty
newdata <- data.frame(factor=1:1)
predict(mod.pdot.u, type='state', newdata=newdata)


# Ditto for the detection probability
coef(mod.pdot.u,    type="det")
SE(mod.pdot.u,      type="det")
confint(mod.pdot.u, type="det")

backTransform(mod.pdot.u, type="det")

predict(mod.pdot.u, type='det', newdata=newdata)

# Look at the posterior probability of detection
# The best unbiased predictor of the posterior probability of detection
bup(ranef(mod.pdot.u))

#----------------------
# p(.) psi(browse)

mod.pdot.psiB.1.u <- unmarked::occu(~1 ~-1+Browsed+Unbrowsed, data=weta.UMF, se=TRUE)
names(mod.pdot.psiB.1.u)

# better to use predict() to do the estimation.
# Because there are no covariates here, the data.frame is empty
newdata <- data.frame(Browsed=c(1,0), Unbrowsed=c(0,1))
cbind(newdata,predict(mod.pdot.psiB.1.u, type='state', newdata=newdata))

coef(mod.pdot.psiB.1.u)
vcov(mod.pdot.psiB.1.u)
mod.pdot.psiB.1.u.oddsratio.browse <- exp( sum(c(1,-1,0)*coef(mod.pdot.psiB.1.u)))
mod.pdot.psiB.1.u.oddsratio.browse

est <- linearComb(mod.pdot.psiB.1.u, c(1,-1), type="state")
exp(est@estimate)


# using the BrowCat variable directly as cell effect models
mod.pdot.psiB.2.u <- unmarked::occu(~1 ~BrowCat, data=weta.UMF, se=TRUE)

# better to use predict() to do the estimation.
# Because there are no covariates here, the data.frame is empty
newdata <- data.frame(BrowCat=c("B","N"))
cbind(newdata,predict(mod.pdot.psiB.2.u, type='state', newdata=newdata))

coef(mod.pdot.psiB.2.u)
vcov(mod.pdot.psiB.2.u)
mod.pdot.psiB.2.u.oddsratio.browse <- exp( sum(c(0,-1,0)*coef(mod.pdot.psiB.2.u)))
mod.pdot.psiB.2.u.oddsratio.browse

est <- linearComb(mod.pdot.psiB.2.u, c(0,-1), type="state")
est
confint(est)
exp(est@estimate)
exp(confint(est))


#---------------------------
# Impact of browse on detectability as well?.

mod.pB.psiB.u <- unmarked::occu(~BrowCat ~BrowCat, data=weta.UMF, se=TRUE)

# better to use predict() to do the estimation.
# Because there are no covariates here, the data.frame is empty
newdata <- data.frame(BrowCat=c("B","N"))
cbind(newdata,predict(mod.pB.psiB.u, type='state', newdata=newdata))

newdata <- data.frame(BrowCat=c("B","N"))
cbind(newdata,predict(mod.pB.psiB.u, type='det', newdata=newdata))

coef(mod.pB.psiB.u)
mod.pB.psiB.u.oddsratio.browse <- exp( sum(c(0,-1,0,0)*coef(mod.pB.psiB.u)))
mod.pB.psiB.u.oddsratio.browse

est <- linearComb(mod.pB.psiB.u, c(0,-1), type="state")
exp(est@estimate)

#-------
# Model where p varies by time psi(browse) 
# 
mod.pt.psiB.u <- unmarked::occu(~Time ~BrowCat, data=weta.UMF, se=TRUE)

# better to use predict() to do the estimation.
# Because there are no covariates here, the data.frame is empty
newdata <- data.frame(BrowCat=c("B","N"))
cbind(newdata,predict(mod.pt.psiB.u, type='state', newdata=newdata))

newdata <- expand.grid(Time=c("T1","T2","T3","T4","T5"))
cbind(newdata,predict(mod.pt.psiB.u, type='det', newdata=newdata))

coef(mod.pt.psiB.u)
mod.pt.psiB.u.oddsratio.browse <- exp( sum(c(0,-1,0,0,0,0,0,0,0)*coef(mod.pt.psiB.u)))
mod.pt.psiB.u.oddsratio.browse

est <- linearComb(mod.pt.psiB.u, c(0,-1), type="state")
exp(est@estimate)

#-------
# Model where p(obs) psi(browse) 
# 
mod.pO.psiB.u <- unmarked::occu(~obs ~BrowCat, data=weta.UMF, se=TRUE)

# better to use predict() to do the estimation.
# Because there are no covariates here, the data.frame is empty
newdata <- data.frame(BrowCat=c("B","N"))
cbind(newdata,predict(mod.pO.psiB.u, type='state', newdata=newdata))

newdata <- expand.grid(obs=c("1",'2','3'))
cbind(newdata,predict(mod.pO.psiB.u, type='det', newdata=newdata))

coef(mod.pO.psiB.u)
mod.pO.psiB.u.oddsratio.browse <- exp( sum(c(0,-1,0,0,0)*coef(mod.pO.psiB.u)))
mod.pO.psiB.u.oddsratio.browse

est <- linearComb(mod.pO.psiB.u, c(0,-1), type="state")
exp(est@estimate)


#-------
# Model where p varies by observer + visit but constant over time psi(browse) 
# 
mod.pOpV.psiB.u <- unmarked::occu(~obs+Time ~BrowCat, data=weta.UMF, se=TRUE)

# better to use predict() to do the estimation.
# Because there are no covariates here, the data.frame is empty
newdata <- data.frame(BrowCat=c("B","N"))
cbind(newdata,predict(mod.pOpV.psiB.u, type='state', newdata=newdata))

newdata <- expand.grid(obs=c("1",'2','3'), Time=c("T1","T2","T3","T4","T5"))
cbind(newdata,predict(mod.pOpV.psiB.u, type='det', newdata=newdata))

coef(mod.pOpV.psiB.u)
mod.pOpV.psiB.u.oddsratio.browse <- exp( sum(c(0,-1,0,0,0,0,0,0,0)*coef(mod.pOpV.psiB.u)))
mod.pOpV.psiB.u.oddsratio.browse

est <- linearComb(mod.pOpV.psiB.u, c(0,-1), type="state")
exp(est@estimate)

#-------
# Other models 

mod.pOpVpB.psiB.u <- unmarked::occu(~obs+Time+BrowCat ~BrowCat, data=weta.UMF, se=TRUE)

mod.pOpVpB.psi.u <- unmarked::occu(~obs+Time+BrowCat ~1, data=weta.UMF, se=TRUE)

mod.pOpVpB.psi.u <-  unmarked::occu(~obs+Time ~1, data=weta.UMF, se=TRUE)


#------
# Model averaging
models.u <-unmarked::fitList(
             mod.pdot.u, 
             mod.pdot.psiB.1.u,
             mod.pB.psiB.u,
             mod.pt.psiB.u,
             mod.pO.psiB.u,
             mod.pOpV.psiB.u,
             mod.pOpVpB.psiB.u,
             mod.pOpVpB.psi.u
             )
aic.table.u <- unmarked::modSel(models.u)
aic.table.u

# Get model averaged estimates of occupancy
predict(models.u, type="state")[1,]

# Get model averaged estimates of detection. 
# Notice these are one big list of nsites x nvisits long
# with the detection probabilities for each site
cbind(survey.covar.u, predict(models.u, type="det"))[1:10,]
