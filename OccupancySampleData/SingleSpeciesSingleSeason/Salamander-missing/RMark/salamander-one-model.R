# Blue Ridge Salamander with some missing values
# Notice how we need to change any NA to . when creating the capture history

# Single Species Single Season Occupancy

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

input.data <- readxl::read_excel("../salamander.xls",
                                 sheet="MissingData",
                                 na="-",
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

# modify the ddl if needed (e.g. for site level covariates)
sal.ddl <- make.design.data(sal.data)
sal.ddl

# What are the parameter names for Single Season Single Species models
setup.parameters("Occupancy", check=TRUE)


# Fit a model
# Note that formula DO NOT HAVE AN = SIGN
mod.fit <-  RMark::mark(sal.data, ddl=sal.ddl,
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

# cleanup
cleanup()

