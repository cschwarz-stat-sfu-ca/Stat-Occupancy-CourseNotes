# Single Species, Multi Season Occupancy analyais

# Northern Spotted Owl (Strix occidentalis caurina) in California.
# s=55 sites visited up to K=5 times per season between 1997 and 2001 (Y=5).
# Detection probabilities relatively constant within years, but likely different among years.

# Using RMark

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

library(readxl)
library(RMark)
library(ggplot2)

# get the data read in
# Data for detections should be a data frame with rows corresponding to sites
# and columns to visits.
# The usual 1=detected; 0=not detected; NA=not visited is used.

input.data <- read.csv(file.path("..","NSO.csv"), header=FALSE, skip=2, na.strings="-")
input.data$V1 <- NULL # drop the site number

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.data)
ncol(input.data)
range(input.data, na.rm=TRUE)
sum(is.na(input.data))


# Extract the history records and create a capture history
input.history <- data.frame(freq=1,
                            ch=apply(input.data,1,paste, collapse=""), stringsAsFactors=FALSE)
head(input.history)


# Change any NA to . in the chapter history
input.history$ch <- gsub("NA",".", input.history$ch, fixed=TRUE)
head(input.history)

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.history)


# Create the data structure
max.visit.per.year <- 8
n.season <- 5

nso.data <- process.data(data=input.history, model="RDOccupEG",
                         time.intervals=c( rep( c(rep(0,max.visit.per.year-1),1),n.season-1),
                                           rep(0,max.visit.per.year-1)))

summary(nso.data)

# What are the parameter names for Single Season Single Species models
setup.parameters("RDOccupEG", check=TRUE)

#there are other parameterizations available
setup.parameters("RDOccupPE", check=TRUE) # psi, epsilon, p
setup.parameters("RDOccupPG", check=TRUE) # psi, gamma, p





#-----
# Fit a models.
# We start with the first formulation in terms of psi, gamma, epsilon, p (do.1_)
# Note that formula DO NOT HAVE AN = SIGN
mod.fit <- RMark::mark(nso.data,
                   model="RDOccupEG",
                   model.parameters=list(
                     Psi   =list(formula=~1), # initial occupancy
                     p     =list(formula=~1),
                     Epsilon=list(formula=~1),
                     Gamma  =list(formula=~1)))

names(mod.fit)
names(mod.fit$results)

# Estimate of initial occupancy
names(mod.fit$results$real)
mod.fit$results$real

get.real(mod.fit, parameter="Psi", se=TRUE)
get.real(mod.fit, parameter="p"  , se=TRUE)
get.real(mod.fit, parameter="Epsilon", se=TRUE)
get.real(mod.fit, parameter="Gamma", se=TRUE)


# Derived parameters - estimated occupancy for each unit in years 2....
names(mod.fit$results$derived)

mod.fit$results$derived$"psi Probability Occupied"
mod.fit$results$derived$"lambda Rate of Change"
mod.fit$results$derived$"log odds lambda"


###################################################
# The random occupancy model is fit as described in
# https://sites.warnercnr.colostate.edu/gwhite/occupancy-estimation-robust-design/
#
# We need to set gamma = 1-epsilon using design matrices
#
#   One model often of interest is the random extinction and colonization model, 
#    where the probability of a site being occupied at time t+1 is the same regardless 
#    of whether or not the site was occupied at time t.  
#   You can obtain estimates for this model in MARK by clever 
#  coding of the design matrix.  

#  Suppose the epsilon of interest is parameter 2 in the PIM, 
#  and the gamma of interest is parameter 3 in the PIM.  
#  Use a common beta parameter for both epsilon and gamma, 
#  i.e., there will be a single column that is modeled by both rows 2 and 3 of the design matrix.  
#   If you specify a logit link for both epsilon and gamma, and code epsilon as -1 
#   in the design matrix and gamma as +1, the result is the model with 1 - epsilon = gamma.  
#   The design matrix looks like:
#    Columns: B1 B2 B3
#  Row 1: 1 0 0
#  Row 2: 0 -1 0 /* This is epsilon */
#  Row 3: 0 1 0 /* This is gamma */
# 
# Unfortunately, RMark does not allow you fit this model directly because only the logit link is
# suported.
#
# You can get the equivalent of a random occupancy model by stacking the data in a special way.
# Contact me for details.



