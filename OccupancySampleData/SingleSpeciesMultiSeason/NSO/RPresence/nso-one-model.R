# Single Species, Multi Season Occupancy analyais

# Northern Spotted Owl (Strix occidentalis caurina) in California.
# s=55 sites visited up to K=5 times per season between 1997 and 2001 (Y=5).
# Detection probabilities relatively constant within years, but likely different among years.

# 2018-08-15 Code contributed by Carl James Schwarz (cschwarz.stat.sfu.cs@gmail.com)

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

#  RPresence package

library(readxl)
library(RPresence)
library(ggplot2)

# Get the RPResence additional functions 
source(file.path("..","..","..","AdditionalFunctions","Rpresence.additional.functions.R"))

# get the data read in
# Data for detections should be a data frame with rows corresponding to sites
# and columns to visits.
# The usual 1=detected; 0=not detected; NA=not visited is used.

input.history <- read.csv(file.path("..","NSO.csv"), header=FALSE, skip=2, na.strings="-")
input.history$V1 <- NULL # drop the site number

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.history)
ncol(input.history)
range(input.history, na.rm=TRUE)
sum(is.na(input.history))


# Five years with 8 visits. Don't need same number of visits/season
Nvisits.per.season  <- rep(8,5) # five years with 8 visits. Don't need same number of visits/season

# Create the *.pao file
nso.pao <- RPresence::createPao(input.history,
                                nsurveyseason=Nvisits.per.season,
                                title='NSO SSMS')
nso.pao

#-----
# Fit a models.
# We start with the first formulation in terms of psi, gamma, epsilon, p (do.1_)
# Note that formula DO NOT HAVE AN = SIGN
mod.fit<- RPresence::occMod(
           model=list(psi~1, gamma~1, epsilon~1, p~1), 
                                        type="do.1", data=nso.pao)
mod.fit <- RPresence.add.derived(mod.fit)  # combine the psi; add the lambda
summary(mod.fit)

names(mod.fit)
names(mod.fit$real)

# Estimate of initial occupancy
mod.fit$real$psi[1,]

# Derived parameters - estimated occupancy for each unit in years 2....
names(mod.fit$derived)
mod.fit$derived$psi[ grepl('unit1_', row.names(mod.fit$derived$psi)),]

# Additional derived parameters - all of the psi stacked together
mod.fit$derived$all_psi[ grepl('unit1_', row.names(mod.fit$derived$all_psi)),]


# Estimate of  local extinction probability for each unit
mod.fit$real$epsilon[ seq(1, by=nrow(input.history), length.out=length(Nvisits.per.season)-1),]

# Estimate of  local colonization probability for each unit
mod.fit$real$gamma[ seq(1, by=nrow(input.history), length.out=length(Nvisits.per.season)-1),]

# Estimate of probability of detection at each time point for each unit
mod.fit$real$p[ grepl('_unit1$', row.names(mod.fit$real$p)),]

# Derived parameters - estimated occupancy for each unit in years 2....
names(mod.fit$derived)
mod.fit$derived$lambda[ grepl('unit1_', row.names(mod.fit$derived$lambda)),]


# Get the change in occupancy
# Not yet possible to estimate the se of these values. May have to use bootstrapping.
mod.fit$derived$lambda [grepl('unit1_', row.names(mod.fit$derived$lambda),  fixed=TRUE),]
mod.fit$derived$lambdap[grepl('unit1_', row.names(mod.fit$derived$lambdap), fixed=TRUE),]



