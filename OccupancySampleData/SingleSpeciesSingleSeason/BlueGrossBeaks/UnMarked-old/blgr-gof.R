# Single Species Single Season Occupancy models 

# Goodness of fit for the gross beak example


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



#---------------------------------
# Do we need to worry about lack of it

# See https://groups.google.com/forum/#!searchin/unmarked/goodness$20of$20fit%7Csort:date/unmarked/3wvnDyLlxok/ct1U723vCAAJ
# https://www.rdocumentation.org/packages/AICcmodavg/versions/2.1-1/topics/mb.gof.test
# MacKenzie, D. I., Bailey, L. L. (2004) Assessing the fit of site-occupancy models. Journal of Agricultural, Biological, and Environmental Statistics 9, 300--318.

library(AICcmodavg)

#  Mackenzie Bailey Goodness of fit test
#  compute chi-square for the top model
mod.pdot.u.chi   = mb.chisq(mod.pdot.u)
mod.pdot.u.chi

mod.pdot.u.boot = mb.gof.test(mod.pdot.u, nsim = 100) # should do about 1000 sims 

print(mod.pdot.u.boot, digit.vals=4, digits.chisq=4)

# estimating chat
names(mod.pdot.u.boot)
chat <- mod.pdot.u.boot$c.hat.est
chat


# Not possible to insert chat in the unmarked::modSel (groan) so we use te aictab from AICmodavg
AICcmodavg::aictab(list(mod.pdot.u, mod.pt.u, mod.pdot.u) )
AICcmodavg::aictab(list(mod.pdot.u, mod.pt.u, mod.pdot.u), c.hat=max(1,chat) )


# The bootstrap goodness of fit can also be done directly in umarked by defining
# our own statistics.
# Function returning three fit-statistics.
fitstats <- function(fm) {
    observed <- getY(fm@data)
    expected <- fitted(fm)
    resids <- residuals(fm)
    sse <- sum(resids^2)
    chisq <- sum((observed - expected)^2 / expected)
    freeTuke <- sum((sqrt(observed) - sqrt(expected))^2)
    out <- c(SSE=sse, Chisq=chisq, freemanTukey=freeTuke)
    return(out)
    }

mod.pdot.u.pboot <- unmarked::parboot(mod.pdot.u, fitstats, nsim=99, report=100)
mod.pdot.u.pboot@t0
head(mod.pdot.u.pboot@t.star)

# How extreme is t0 relative to the bootstrap values.
# Small probabilities indicate lack of fit
colMeans(mod.pdot.u.pboot@t.star > matrix(mod.pdot.u.pboot@t0, byrow=TRUE,
                                    ncol=length(mod.pdot.u.pboot@t0),
                                    nrow=nrow(mod.pdot.u.pboot@t.star)))

# No easy way to automatically inflate the se and ci using c.hat in the unmarked package
