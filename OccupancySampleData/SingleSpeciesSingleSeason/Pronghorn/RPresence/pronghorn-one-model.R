# Example using model averaging on the pronghorn dataset
# 256 sites, two sampling occasions per site
# Covariates used in this example:
# 1. sagebrush (continuos) - Sagebrush density
# 2. aspect (Categorical) - Compass direction slope faces (N,S,E,W)

# Code contributed by Neil Faught - 26/11/2018

#  RPresence package

library(readxl)
library(RPresence)
library(ggplot2)

# get the data read in
# Data for detections should be a data frame with rows corresponding to sites
# and columns to visits.
# The usual 1=detected; 0=not detected; NA=not visited is used.

input.data<-read.csv(file.path("..","pronghorn.csv"), header=TRUE, as.is=TRUE, strip.white=TRUE)
head(input.data)

input.history <- input.data[,c(2,3)] # the history extracted
head(input.history)

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.history)
ncol(input.history)
range(input.history, na.rm=TRUE) # check that all values are either 0 or 1
sum(is.na(input.history))    # are there any missing values?

# Get the sagebrush and aspect information. These are the covariates that will
# be used in the analysis. These are unit-level covariates.
unit.cov = input.data[,c(4,7)]
head(unit.cov)

# Create the *.pao file
prong.pao <- RPresence::createPao(input.history,
                                  unitcov=unit.cov,
                                  title='Pronghorn SSSS')
prong.pao

# Model where occupancy varies with sagebrush density, and detection probability
# varies with aspect
mod.fit <- RPresence::occMod(model=list(psi~sagebrush, p~aspect), 
                              type="so", data=prong.pao)
summary(mod.fit)
head(mod.fit$real$psi)
mod.fit$real$p[seq(1, by=nrow(input.history), length.out=ncol(input.history)),]

# plot occupancy as a function of elevation of stream
plotdata = data.frame(mod.fit$real$psi, sagebrush = unit.cov$sagebrush)

ggplot(data=plotdata, aes(x=sagebrush, y=est))+
  ggtitle("Occupancy as a function of sagebrush density")+
  geom_point()+
  geom_ribbon(aes(ymin=lower_0.95, ymax=upper_0.95), alpha=0.2)+
  ylim(0,1)+
  ylab("Estimated occupancy")


