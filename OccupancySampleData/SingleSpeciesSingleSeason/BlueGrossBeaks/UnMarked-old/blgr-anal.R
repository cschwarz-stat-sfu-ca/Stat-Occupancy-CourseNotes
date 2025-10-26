# Single Species Single Season Occupancy 

# Blue Gross Beaks.
#Downloaded from https://sites.google.com/site/asrworkshop/home/schedule/r-occupancy-1

#An occupancy study was made on Blue Grosbeaks (Guiraca caerulea) 
# on 41 old fields planted to longleaf pines (Pinus palustris) 
# in southern Georgia, USA. 

# Surveys were 500 m transects across each field 
# and were completed three times during the breeding season in 2001.

# Columns in the file are:
#    field - field number
#    v1, v2, v3 -  detection histories for each site on each of 3 visit during the 2001 breeding season.    
#    field.size - size of the files
#    bqi - Enrollment in bobwihte quail initiative; does occupancy increase if field belongs to this initiative?
#    crop.hist - crop history
#    crop1, crop2 - indicator variables for the crop history
#    count1, count2, count3 - are actual counts of birds detected in each visit

# 2018-08-15 Code contributed by Carl James Schwarz (cschwarz.stat.sfu.cs@gmail.com)




#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

#  unmarked package

library(readxl)
library(unmarked)
library(ggplot2)

# get the data read in
# Data for detections should be a data frame with rows corresponding to sites
# and columns to visits.
# The usual 1=detected; 0=not detected; NA=not visited is used.

input.data <- read.csv(file.path("..","blgr.csv"), 
                  header=TRUE, as.is=TRUE, strip.white=TRUE) 
head(input.data)

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.data)
range(input.data[, c("v1","v2","v3")], na.rm=TRUE)
sum(is.na(input.data[, c("v1","v2","v3")]))

input.history <- input.data[, c("v1","v2","v3")]
head(input.history)


# Create the data structore for the occupancy model
# unmarked does not have any reserved keywords like RPresence so 
# we need to construct a covariate data frame with the covariates
# again stacked. The ordering is different in unmarked.
# Here you want for site 1, the visit specific covariates, then
# for site 2, etc.
# This needs to be added to the unmarked dataframe.

# We also want to set up covariate values where we
# set the detection probability equal in first 2 occasions and last 3 occasions.
# We need to define a survey covariate that has two levels 
obs.covar <- data.frame( Site=rep(1:nrow(input.history), each=ncol(input.history)),
                         Visit=rep(1:ncol(input.history), nrow(input.history)),
                         Time =rep(c("T1","T2","T3"), nrow(input.history)))
obs.covar[1:10,]

grossbeak.UMF <- unmarked::unmarkedFrameOccu(
                     y = input.history,
                     obsCovs=obs.covar) 
summary(grossbeak.UMF)
plot(grossbeak.UMF)



#-----
# Fit some models.
# Note that formula DO NOT HAVE AN = SIGN
# The two formula are for detecton and occupancy in that order
mod.pdot.u <- unmarked::occu(~1 ~1 , data=grossbeak.UMF, se=TRUE)

# This creates a complicated S4 object with SLOTS that can be accessed using various
# functions. This is more complcated than with RPresence
slotNames(mod.pdot.u)

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


#-------
# Model where p(t) varies across survey occasions
# 

# Note that formula DO NOT HAVE AN = SIGN
# The two formula are for detecton and occupancy in that order
mod.pt.u <- unmarked::occu(~Time ~1 , data=grossbeak.UMF, se=TRUE)

# This creates a complicated S4 object with SLOTS that can be accessed using various
# functions. This is more complcated than with RPresence

# These functions return estimates on the LOGIT scale - not very useful for most people
# better to use predict() to do the estimation.
# Because there are no covariates here, the data.frame is empty
newdata <- data.frame(Time=c("T1","T2","T3"), stringsAsFactors=FALSE)
predict(mod.pt.u, type='state', newdata=newdata)[1,]

# Ditto for the detection probability
cbind(newdata, predict(mod.pt.u, type='det', newdata=newdata))

# Look at the posterior probability of detection
# The best unbiased predictor of the posterior probability of detection
bup(ranef(mod.pt.u))
sum(bup(ranef(mod.pt.u)))


#------
# Model averaging
models.u <-unmarked::fitList(mod.pdot.u, mod.pt.u)
aic.table.u <- unmarked::modSel(models.u)
aic.table.u

# Get model averaged estimates of occupancy
predict(models.u, type="state")[1,]

# Get model averaged estimates of detection. 
# Notice these are one big list of nsites x nvisits long
# with the detection probabilities for each site
cbind(obs.covar, predict(models.u, type="det"))[1:10,]

