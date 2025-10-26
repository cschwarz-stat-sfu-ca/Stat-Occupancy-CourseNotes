# Single Species Single Season Occupancy 

# Salamander Example with model averaging with individual model fits
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

input.data <- readxl::read_excel("../salamander.xls",
                                 sheet="CompleteData",
                                 col_names=FALSE)  # notice no column names in row 1 of data file. 

head(input.data)


# Extract the history records and create a capture history
input.history <- data.frame(freq=1,
                            ch=apply(input.data[,1:5],1,paste, collapse=""), stringsAsFactors=FALSE)
head(input.history)

# Change any NA to . in the chapter history
input.history$ch <- gsub("NA",".", input.history$ch, fixed=TRUE)
head(input.history)

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.history)



sal.data <- process.data(data=input.history,
                          model="Occupancy")

summary(sal.data)

# Fit some models.
# Notice that RMark does not allow missing values in time-varying covariates, even when visits are not made
mod.fit <-  RMark::mark(sal.data,
                        model="Occupancy",
                        model.parameters=list(
                          Psi   =list(formula=~1),
                          p     =list(formula=~1) 
                        )
)
summary(mod.fit)

# Look the objects returned in more details
names(mod.fit)

# look at estimates on beta and original scale
mod.fit$results$beta  # on the logit scale

mod.fit$results$real# on the regular 0-1 scale for each site

# derived variables is the occupancy probability 
names(mod.fit$results$derived)
mod.fit$results$derived$Occupancy

# Model where p(t) varies across survey occasions
# 
mod.fit.pt <-  RMark::mark(sal.data,
                        model="Occupancy",
                        model.parameters=list(
                          Psi   =list(formula=~1),
                          p     =list(formula=~time) 
                        )
)
summary(mod.fit.pt)

# Look the objects returned in more details
names(mod.fit.pt)

# look at estimates on beta and original scale
mod.fit.pt$results$beta  # on the logit scale

mod.fit.pt$results$real# on the regular 0-1 scale for each site

# derived variables is the occupancy probability 
names(mod.fit.pt$results$derived)
mod.fit.pt$results$derived$Occupancy

# Fit a model with the detection probability equal in first 2 occasions and last 3 occasions.
# We need to define a survey covariate that has two levels 
# add a survey level covariate. by adding a column to the design matrix
sal.ddl <- make.design.data(sal.data)
sal.ddl

sal.ddl$p$Effort <- factor(c("d1","d1","d2","d2","d2"))
sal.ddl

mod.fit.pcustom <-  RMark::mark(sal.data, ddl = sal.ddl,
                                model="Occupancy",
                                model.parameters=list(
                                  Psi   =list(formula=~1),
                                  p     =list(formula=~Effort) 
                                )
)
summary(mod.fit.pcustom)

# Look the objects returned in more details
names(mod.fit.pcustom)

# look at estimates on beta and original scale
mod.fit.pcustom$results$beta  # on the logit scale

mod.fit.pcustom$results$real# on the regular 0-1 scale for each site

# derived variables is the occupancy probability 
names(mod.fit.pcustom$results$derived)
mod.fit.pcustom$results$derived$Occupancy

# Model averaging
model.set <- RMark::collect.models( type="Occupancy")
model.set

Psi.ma <- RMark::model.average(model.set, param="Psi")
Psi.ma

# cleanup
cleanup(ask=FALSE)

