# Example using model averaging on the pronghorn dataset
# 256 sites, two sampling occasions per site
# Covariates used in this example:
# 1. sagebrush (continuos) - Sagebrush density
# 2. aspect (Categorical) - Compass direction slope faces (N,S,E,W)

# Code contributed by Neil Faught - 26/11/2018

#  RMark package

library(readxl)
library(RMark)
library(ggplot2)

# get the data read in
# Data for detections should be a data frame with rows corresponding to sites
# and columns to visits.
# The usual 1=detected; 0=not detected; NA=not visited is used.

input.data<-read.csv(file.path("..","pronghorn.csv"), header=TRUE, as.is=TRUE, strip.white=TRUE)
head(input.data)

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.data)
range(input.data[, c("Survey.1","Survey.2")], na.rm=TRUE)
sum(is.na(input.data[, c("Survey.1","Survey.2")]))

input.history <- input.data[, c("Survey.1","Survey.2")]
head(input.history)

site.covar <- input.data[, c("sagebrush","aspect")]
head(site.covar)
#Convert aspect to a factor
site.covar$aspect = as.factor(site.covar$aspect)

# Extract the history records and create a capture history
input.history <- data.frame(freq=1,
                            ch=apply(input.history,1,paste, collapse=""), 
                            stringsAsFactors=FALSE)
head(input.history)

# Change any NA to . in the capture history
input.history$ch <- gsub("NA",".", input.history$ch, fixed=TRUE)
head(input.history)

#Add aspect and sagebrush data to the capture history
input.history = cbind(input.history, site.covar)
head(input.history)

#Create the data structure
pronghorn.data <- process.data(data = input.history,
                               group="aspect",
                               model = "Occupancy")
summary(pronghorn.data)

# modify the ddl if needed (e.g. for survey/visit level covariates)
pronghorn.ddl <- make.design.data(pronghorn.data)
pronghorn.ddl

# What are the parameter names for Single Season Single Species models
setup.parameters("Occupancy", check=TRUE)

# Fit a model
# Model where detection probability (p) varies with aspect, and
# Occupancy (Psi) varies with sagebrush
# Notice that RMark does not allow missing values in time-varying covariates, even when visits are not made

mod.fit <-  RMark::mark(pronghorn.data, ddl = pronghorn.ddl,
                        model="Occupancy",
                        model.parameters=list(
                          Psi   =list(formula=~sagebrush),
                          p     =list(formula=~aspect) 
                        )
)

summary(mod.fit)

# Look at the objects returned in more details
names(mod.fit)

# look at estimates on beta and original scale
mod.fit$results$beta  # on the logit scale

mod.fit$results$real# on the regular 0-1 scale for each site

# derived variables is the occupancy probability 
names(mod.fit$results$derived)

mod.fit$results$derived$Occupancy

# get the psi-values as function of the covariates
pronghorn.ddl$Psi # see the index numbers

sage.df <-data.frame(sagebrush=seq(min(input.history$sagebrush,na.rm=TRUE),
                                   max(input.history$sagebrush, na.rm=TRUE),.5)) 
pred.psi <- covariate.predictions(mod.fit, indices=9, data=sage.df)
head(pred.psi$estimates)

ggplot(data=pred.psi$estimates, aes(x=covdata, y=estimate))+
  ggtitle("Occupancy probability as a function of sagebrush density")+
  geom_point()+
  geom_ribbon(aes(ymin=lcl, ymax=ucl), alpha=.2)+
  ylim(0,1)+
  ylab("Estimated probability of occupancy")+
  xlab("Sagebrush")

# It is often convenient to estimate occupancy for each site using the
# values of the covariates specific for that site.
# This is similar to RPresence which gives estimates at each site.

# We create a data frame with the covariates etc as measure on the input.history dataframe.
# But we also need to add the appropriate index number for psi to the covariate data frame to account
# for the groups definition
pronghorn.ddl$Psi
site.covar <- merge(pronghorn.ddl$Psi, site.covar)
head(site.covar)

# we only want a prediction for that particular model.index.
# we need to create a variable "index" that has the model.index
site.covar$index <- site.covar$model.index

site.pred <- covariate.predictions(mod.fit,
                                   data=site.covar)
ggplot(site.pred$estimates, aes(x=sagebrush, y=estimate))+
  ggtitle("Model predictions of occupancy for each site")+
  geom_point()+
  geom_ribbon(aes(ymin=lcl, ymax=ucl), alpha=0.2)

#cleanup
cleanup(ask=FALSE)
