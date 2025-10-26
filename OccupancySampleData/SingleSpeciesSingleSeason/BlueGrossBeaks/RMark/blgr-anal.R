# Single Species Single Season Occupancy models using RMark 

# Blue Gross Beaks.
#Downloaded from https://sites.google.com/site/asrworkshop/home/schedule/r-occupancy-1

#An occupancy study was made on Blue Grosbeaks (Guiraca caerulea) 
# on 41 old fields planted to longleaf pines (Pinus palustris) 
# in southern Georgia, USA. 

# Surveys were 500 m transects across each field 
# and were completed three times during the breeding season in 2001.

# Columns in the file are:
#    field - field number
#    v1, v2, v3 -  detection histories for each site on each of 3 visit during the 2001 breeding season.    
#    field.size - size of the files
#    bqi - Enrollment in bobwihte quail initiative; does occupancy increase if field belongs to this initiative?
#    crop.hist - crop history
#    crop1, crop2 - indicator variables for the crop history
#    count1, count2, count3 - are actual counts of birds detected in each visit


# 2018-08-15 Code contributed by Carl James Schwarz (cschwarz.stat.sfu.cs@gmail.com)


#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

#  RMark package
library(readxl)
library(RMark)
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

# Extract the history records and create a capture history
input.history <- data.frame(freq=1,
                            ch=apply(input.data[,c("v1","v2","v3")],1,paste, collapse=""), 
                            stringsAsFactors=FALSE)
head(input.history)

# Change any NA to . in the chapter history
input.history$ch <- gsub("NA",".", input.history$ch, fixed=TRUE)
head(input.history)

#Add remaining columns (field.size, crop.hist, etc) to capture history
#These rows are contained in the last 10 columns
input.history = cbind(input.history, input.data[,5:14])
#Remove empty column
input.history$X <- NULL

#Create the data structure
grossbeak.data <- process.data(data = input.history,
                               model = "Occupancy")
summary(grossbeak.data)

# If survey covariates are present, modify the ddl
grossbeak.ddl <- make.design.data(grossbeak.data)
grossbeak.ddl


# What are the parameter names for Single Season Single Species models
setup.parameters("Occupancy", check=TRUE)

# Fit a model
# Notice that RMark does not allow missing values in time-varying covariates, even when visits are not made
mod.fit <-  RMark::mark(grossbeak.data, ddl = grossbeak.ddl,
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

#-------
# Model where p(t) varies across survey occasions
# 

mod.fit.pt <-  RMark::mark(grossbeak.data, ddl = grossbeak.ddl,
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


#------
# Model averaging

# collect models and make AICc table

model.set <- RMark::collect.models( type="Occupancy")
model.set

names(model.set)
model.set$model.table

# model averaged values
Psi.ma <- RMark::model.average(model.set, param="Psi")
Psi.ma

p.ma <- RMark::model.average(model.set, param="p")
p.ma

#cleanup
cleanup(ask=FALSE)