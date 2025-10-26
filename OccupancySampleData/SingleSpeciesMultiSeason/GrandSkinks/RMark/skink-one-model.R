# Single Species, MultiSeason Occupancy analyais

# Grand Skinks

#   Data has been collected on 352 tors over a 5 year (the seasons) period,^
#   although not all tors (rock piles) were surveyed each year, with up to 3 surveys^
#   of each tor per year.  

#   The 15 columns are in 5 blocks of 3.

#   There is also a site-specific covariate Pasture indicating whether the surrounding^
#   matrix is either predominately the modified habitat (farm pasture, Pasture =1) or
#   "native" grassland (tussock, Pasture = 0).

# 2018-11-09 Code submitted by Neil Faught

# Using the RMark package

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

library(car)
library(readxl)
library(RMark)
library(ggplot2)

# Get the RMark additional functions 
source(file.path("..","..","..","AdditionalFunctions","RMark.additional.functions.R"))

# get the data read in
# Data for detections should be a data frame with rows corresponding to sites
# and columns to visits.
# The usual 1=detected; 0=not detected; NA=not visited is used.

input.data <- readxl::read_excel(file.path("..","GrandSkinks.xls"), 
                                 sheet="DetectionHistory",
                                 skip=3, na="-",
                                 col_names=FALSE)
head(input.data)

input.history <- input.data[, 2:16]
head(input.history)

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.history)
ncol(input.history)
range(input.history, na.rm=TRUE)
#Looks like there are some missing values
sum(is.na(input.history))

site.covar <- data.frame(Site=1:nrow(input.data),
                         Pasture=car::recode(input.data$X__18,
                                             "1='Modified'; 0='Native';"))
head(site.covar)

#Format the capture history to be used by RMark
input.history <- data.frame(freq=1,
                            ch=apply(input.history,1,paste, collapse=""), stringsAsFactors=FALSE)
head(input.history)

# Change any NA to . in the chapter history
input.history$ch <- gsub("NA",".", input.history$ch, fixed=TRUE)
head(input.history)

#Add site covariates to input history
input.history = cbind(input.history,site.covar)
head(input.history)

#Create the RMark data structure
#Five  Seasons, with 3 visits per season
max.visit.per.year <- 3
n.season <- 5

skink.data <- process.data(data=input.history, 
                           model="RDOccupEG",
                           time.intervals=c( rep( c(rep(0,max.visit.per.year-1),1),n.season-1),
                                             rep(0,max.visit.per.year-1)))
summary(skink.data)

# add visit level covariates
# In this case none
skink.ddl <- make.design.data(skink.data)
skink.ddl

# Fit some models.
# We start with the first formulation in terms of psi, gamma, epsilon, p (RDOccupEG)
# Note that formula DO NOT HAVE AN = SIGN
mod.fit <-  RMark::mark(skink.data, ddl = skink.ddl,
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
# estimated ratio of in occupancy from year to year, and estimated log
# of the ratio of odds in occupancy from year to year
mod.fit$results$derived

mod.fit$results$derived$"psi Probability Occupied"
mod.fit$results$derived$"lambda Rate of Change"
mod.fit$results$derived$"log odds lambda"

# Plot derived occupancy probabilities
plotdata = mod.fit$results$derived$"psi Probability Occupied"
plotdata$Season = c(1:5)
head(plotdata)


ggplot(data=plotdata, aes(x=Season,y=estimate))+
  ggtitle("Estimated occupancy over time")+
  geom_point(position=position_dodge(w=0.2))+
  geom_line(position=position_dodge(w=0.2))+
  ylim(0,1)+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.1,position=position_dodge(w=0.2))+
  scale_x_continuous(breaks=1:10)

cleanup(ask=FALSE)
