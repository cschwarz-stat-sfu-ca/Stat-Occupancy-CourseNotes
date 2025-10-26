

# Single Species Single Season Occupancy models 

# MacKenzie-Bailey goodness of fit test
# Illustration of using bootstrapping to estimate standard errors rather than using C-hat

# Using the unmarked package


# 2018-08-15 Code contributed by Carl James Schwarz (cschwarz.stat.sfu.cs@gmail.com)

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

#  unmarked package

library(AICcmodavg)
library(readxl)
library(unmarked)
library(ggplot2)

# get the data read in
# Data for detections should be a data frame with rows corresponding to sites
# and columns to visits.
# The usual 1=detected; 0=not detected; NA=not visited is used.

input.data <- readxl::read_excel("../salamander.xls",
                                 sheet="CompleteData",
                                 na="-",
                                 col_names=FALSE)  # notice no column names in row 1 of data file. 
head(input.data)

input.history <- input.data[, 1:5]
# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.history)
ncol(input.history)
range(input.history, na.rm=TRUE)
sum(is.na(input.history))

head(input.history)


# Create the data structore for the occupancy model
# unmarked does not have any reserved keywords like RPresence so 
# we need to construct a covariate data frame with the covariates
# again stacked. The ordering is different in unmarked than in Rpresence or Presence.
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
# Fit one model.
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

#-------
# Model where p(t) varies across survey occasions
# 

# Note that formula DO NOT HAVE AN = SIGN
# The two formula are for detecton and occupancy in that order
mod.pt.u <- unmarked::occu(~Time ~1 , data=salamander.UMF, se=TRUE)

# This creates a complicated S4 object with SLOTS that can be accessed using various
# functions. This is more complcated than with RPresence

# These functions return estimates on the LOGIT scale - not very useful for most people
# better to use predict() to do the estimation.
# Because there are no covariates here, the data.frame is empty
newdata <- data.frame(Time=c("T1","T2","T3","T4","T5"))
predict(mod.pt.u, type='state', newdata=newdata)[1,]

# Ditto for the detection probability
cbind(newdata, predict(mod.pt.u, type='det', newdata=newdata))

# Look at the posterior probability of detection
# The best unbiased predictor of the posterior probability of detection
bup(ranef(mod.pt.u))
sum(bup(ranef(mod.pt.u)))


#-----


# Note that formula DO NOT HAVE AN = SIGN
# The two formula are for detecton and occupancy in that order
mod.pcustom.u <- unmarked::occu(~D ~1 , data=salamander.UMF, se=TRUE)

# This creates a complicated S4 object with SLOTS that can be accessed using various
# functions. This is more complcated than with RPresence

# These functions return estimates on the LOGIT scale - not very useful for most people
# better to use predict() to do the estimation.
# Because there are no covariates here, the data.frame is empty
newdata <- data.frame(D=c("D1","D1","D2","D2","D2"))
predict(mod.pcustom.u, type='state', newdata=newdata)[1,]

# Ditto for the detection probability
cbind(newdata, predict(mod.pcustom.u, type='det', newdata=newdata))

# Look at the posterior probability of detection
# The best unbiased predictor of the posterior probability of detection
bup(ranef(mod.pcustom.u))
sum(bup(ranef(mod.pcustom.u)))


#------
# Model averaging
models.u <-unmarked::fitList(mod.pdot.u, mod.pt.u, mod.pcustom.u)
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
mod.pcustom.u.chi   = mb.chisq(mod.pcustom.u)
mod.pcustom.u.chi

mod.pcustom.u.boot = mb.gof.test(mod.pcustom.u, nsim = 100) # should do about 1000 sims 

print(mod.pcustom.u.boot, digit.vals=4, digits.chisq=4)

# estimating chat
names(mod.pcustom.u.boot)
chat <- mod.pcustom.u.boot$c.hat.est
chat


# Not possible to insert chat in the unmarked::modSel (groan) so we use te aictab from AICmodavg
AICcmodavg::aictab(list(mod.pdot.u, mod.pt.u, mod.pcustom.u) )
AICcmodavg::aictab(list(mod.pdot.u, mod.pt.u, mod.pcustom.u), c.hat=chat )


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

mod.pcustom.u.pboot <- unmarked::parboot(mod.pcustom.u, fitstats, nsim=99, report=100)
mod.pcustom.u.pboot@t0
head(mod.pcustom.u.pboot@t.star)

# How extreme is t0 relative to the bootstrap values.
# Small probabilities indicate lack of fit
colMeans(mod.pcustom.u.pboot@t.star > matrix(mod.pcustom.u.pboot@t0, byrow=TRUE,
                                    ncol=length(mod.pcustom.u.pboot@t0),
                                    nrow=nrow(mod.pcustom.u.pboot@t.star)))

# No easy way to automatically inflate the se and ci using c.hat in the unmarked package
