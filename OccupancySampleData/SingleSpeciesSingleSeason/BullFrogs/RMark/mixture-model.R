# Bullfrogs
# Single Species Single Season Occupancy - Mixture model

# Fitting a single model

# 2018-08-15 Code contributed by Carl James Schwarz (cschwarz.stat.sfu.cs@gmail.com)

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

input.data <- read.table(file.path("..","MARK","bullfrog.inp"), as.is=TRUE, colClasses="character")
head(input.data)


# Extract the history records and create a capture history
input.history <- data.frame(freq=1,
                            ch=input.data$V1, stringsAsFactors=FALSE)
head(input.history)

# Change any NA to . in the chapter history
input.history$ch <- gsub("NA",".", input.history$ch, fixed=TRUE)
head(input.history)

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.history)


bullfrog.data <- process.data(data=input.history,
                          model="Occupancy")
summary(bullfrog.data)

# What are the parameter names for Single Season Single Species models
setup.parameters("Occupamcut", check=TRUE)


# Fit a model without heterogeneity
mod.fit <-  RMark::mark(bullfrog.data,
                            model="Occupancy",
                            model.parameters=list(
                              Psi   =list(formula=~1),
                              p     =list(formula=~1)
                             )
)
summary(mod.fit)

# Look the objects returned in more details
names(mod.fit)
names(mod.fit$results)

# look at estimates on beta and original scale
mod.fit$results$beta  # on the logit scale

mod.fit$results$real# on the regular 0-1 scale for each site

# derived variabldes is the occupancy probability 
names(mod.fit$results$derived)

mod.fit$results$derived$Occupancy


# alternatively
get.real(mod.fit, "Psi", se=TRUE)
get.real(mod.fit, "Psi", pim=TRUE)


get.real(mod.fit, "p", se=TRUE)



# Fit a model allowing for heterogeneity
# Note that formula HAVE AN = SIGN

bullfrog.data <- process.data(data=input.history,
                              model="OccupHet")
summary(bullfrog.data)

# What are the parameter names for Single Season Single Species models
setup.parameters("OccupHet", check=TRUE)


mod.fit.het <-  RMark::mark(bullfrog.data,
                        model="OccupHety",   # notice change of modelnames
                        model.parameters=list(
                          Psi   =list(formula=~1),
                          p     =list(formula=~mixture),
                          pi    =list(formula=~1)
                        )
                     )
summary(mod.fit.het)

# Look the objects returned in more details
names(mod.fit.het)
names(mod.fit.het$results)

# look at estimates on beta and original scale
mod.fit.het$results$beta  # on the logit scale

mod.fit.het$results$real# on the regular 0-1 scale for each site

# derived variabldes is the occupancy probability 
names(mod.fit.het$results$derived)

mod.fit.het$results$derived$Occupancy


# alternatively
get.real(mod.fit.het, "Psi", se=TRUE)
get.real(mod.fit.het, "Psi", pim=TRUE)


get.real(mod.fit.het, "p", se=TRUE)



# Fit a Royle Nichols Poisson model

bullfrog.data <- process.data(data=input.history,
                              model="OccupRNPoisson")
summary(bullfrog.data)

# What are the parameter names for Single Season Single Species models
setup.parameters("OccupRNPoisson", check=TRUE)


mod.fit.rp <-  RMark::mark(bullfrog.data,
                            model="OccupRNPoisson",   # notice change of modelnames
                            model.parameters=list(
                              r  =list(formula=~1),
                              Lambda     =list(formula=~1))
                            )

summary(mod.fit.rp)

# Look the objects returned in more details
names(mod.fit.rp)
names(mod.fit.rp$results)

# look at estimates on beta and original scale
mod.fit.rp$results$beta  # on the logit scale

mod.fit.rp$results$real# on the regular 0-1 scale for each site

# derived variabldes is the occupancy probability 
names(mod.fit.rp$results$derived)

mod.fit.rp$results$derived$Occupancy




collect.models()

# cleanup
cleanup(ask=FALSE)

