# Single Species Single Season Occupancy models - with missing data

# Salamander data

# Using unmarked packages


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

input.data <- readxl::read_excel("../salamander.xls",
                                 sheet="MissingData",
                                 na="-",
                                 col_names=FALSE)  # notice no column names in row 1 of data file. 

head(input.data)


# Extract the history records
input.history <- input.data[, 1:5] # the history extracted
head(input.history)

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.history)
ncol(input.history)
range(input.history, na.rm=TRUE) # check that all values are either 0 or 1
sum(is.na(input.history))    # are there any missing values?


# Create the data structore for the occupancy model
# unmarked does not have any reserved keywords like RPresence so 
# we need to construct a covariate data frame with the covariates
# again stacked. The ordering in unmarked.
# Here you want for site 1, the visit specific covariates, then
# for site 2, etc.
# This needs to be added to the unmarked dataframe.

# We also want to set up covariate values where we
# set the detection probability equal in first 2 occasions and last 3 occasions.
# We need to define a survey covariate that has two levels 
obs.covar <- data.frame( Site=rep(1:nrow(input.history), each=ncol(input.history)),
                         Visit=rep(1:ncol(input.history), nrow(input.history)),
                         Time =rep(c("T1","T2","T3","T4","T5"), nrow(input.history)),
                         D =rep(c("D1","D1","D2","D2","D2"), nrow(input.history)))
obs.covar[1:10,]

salamander.UMF <- unmarked::unmarkedFrameOccu(
                     y = input.history,
                     obsCovs=obs.covar) 
summary(salamander.UMF)
plot(salamander.UMF)



#-----
# Fit a model
# Note that formula DO NOT HAVE AN = SIGN
# The two formula are for detecton and occupancy in that order
mod.pdot.u <- unmarked::occu(~1 ~1 , data=salamander.UMF, se=TRUE)

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


