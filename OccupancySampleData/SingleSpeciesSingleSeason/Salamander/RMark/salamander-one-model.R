# Blue Ridge Salamander 
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
                                 sheet="CompleteData",
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

# What are the parameter names for Single Season Single Species models
setup.parameters("Occupancy", check=TRUE)


# Fit a model
# Note that formula HAVE AN = SIGN
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



#######################################
# Fit a model with p varying over visits

mod.fit2 <-  RMark::mark(sal.data,
                        model="Occupancy",
                        model.parameters=list(
                          Psi   =list(formula=~1),
                          p     =list(formula=~time)
                        )
)
summary(mod.fit2)

get.real(mod.fit2, "p", se=TRUE)

####################################################
# fit a model with p equal in first two visits and last 3 visits
# This is a survey-level covariate and so you need to modify the ddl
#

# add a survey level covariate. by adding a column to the design matrix

sal.ddl <- make.design.data(sal.data)
sal.ddl

sal.ddl$p$Effort <- factor(c("d1","d1","d2","d2","d2"))
sal.ddl

mod.fit3 <-  RMark::mark(sal.data, ddl=sal.ddl,
                         model="Occupancy",
                         model.parameters=list(
                           Psi   =list(formula=~1),
                           p     =list(formula=~Effort)
                         )
)
summary(mod.fit3)

get.real(mod.fit3, "p", se=TRUE)


##############################################
# Collect models
model.set <- RMark::collect.models( type="Occupancy")
model.set

names(model.set)
model.set$model.table


# model averaged values
get.real(mod.fit , "Psi", se=TRUE)
get.real(mod.fit2, "Psi", se=TRUE)
get.real(mod.fit3, "Psi", se=TRUE)

Psi.ma <- RMark::model.average(model.set, param="Psi")
Psi.ma


# cleanup
cleanup(ask=FALSE)

