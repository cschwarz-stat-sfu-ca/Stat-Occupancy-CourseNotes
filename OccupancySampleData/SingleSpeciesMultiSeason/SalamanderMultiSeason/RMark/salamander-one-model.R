# Single Species, Multi Season Occupancy analyais

#   Salamanders over four seasons with covariates on elevation and stream type.
#   We will use only the last 2 years (unsure why)

# 2018-08-15 Code contributed by Carl James Schwarz (cschwarz.stat.sfu.cs@gmail.com)

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

#  RPresence package

library(readxl)
library(RMark)
library(ggplot2)

# Get the RMark additional functions 
source(file.path("..","..","..","AdditionalFunctions","RMark.additional.functions.R"))

# get the data read in
# Data for detections should be a data frame with rows corresponding to sites
# and columns to visits.
# The usual 1=detected; 0=not detected; NA=not visited is used.

input.data <- readxl::read_excel(file.path("..","SalamanderMultiSeason.xls"), 
                            sheet="Detections",
                            skip=2, na="-",
                            col_names=FALSE)
head(input.data)



input.history <- input.data[, 10:19]
head(input.history)

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.history)
ncol(input.history)
range(input.history, na.rm=TRUE)
sum(is.na(input.history))

site.covar <- data.frame(Site=1:nrow(input.data),
                         Elevation=unlist(input.data[,21]),
                         Prox=car::recode(unlist(input.data[,22]),
                                            "1='Yes'; 0='No';"))
row.names(site.covar) <- NULL
head(site.covar)

#Format the capture history to be used by RMark
input.history <- data.frame(freq=1,
                            ch=apply(input.history,1,paste, collapse=""), stringsAsFactors=FALSE)
head(input.history)

#Add site covariates to input history
input.history = cbind(input.history,site.covar)
head(input.history)

#Create the RMark data structure
#Two  Seasons, with 5 visits per season
max.visit.per.year <- 5
n.season <- 2

salamander.data <- process.data(data=input.history, 
                           model="RDOccupEG",
                           time.intervals=c( rep( c(rep(0,max.visit.per.year-1),1),n.season-1),
                                             rep(0,max.visit.per.year-1)))
summary(salamander.data)

# Fit a models.
# We start with the first formulation in terms of psi, gamma, epsilon, p
mod.fit <-  RMark::mark(salamander.data,
                        model="RDOccupEG",
                        model.parameters=list(
                          Psi   =list(formula=~1),
                          p     =list(formula=~1),
                          Epsilon = list(formula=~1),
                          Gamma = list(formula=~1)
                        )
)
summary(mod.fit)
# Estimates of real initial occupancy, detection, local exctinction,
# and local colonization probabilities
mod.fit$results$real

# Estimates of regression coefficients
mod.fit$results$beta

# Estimates of various derived parameters for estimated occpancy in each year,
# estimated ratio of in occupancy from year to year, and estimated log of the
# of the ratio of odds of  occupancy from year to year
mod.fit$results$derived

mod.fit$results$derived$"psi Probability Occupied"
mod.fit$results$derived$"lambda Rate of Change"
mod.fit$results$derived$"log odds lambda"

cleanup(ask = FALSE)
