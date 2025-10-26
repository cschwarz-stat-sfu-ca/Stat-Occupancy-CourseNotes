# Single Species, MultiSeason Occupancy analyais

# Grand Skinks

#   Ships with PRESENCE.

#   Data has been collected on 352 tors over a 5 year (the seasons) period,^
#   although not all tors (rock piles) were surveyed each year, with up to 3 surveys^
#   of each tor per year.  

#   The 15 columns are in 5 blocks of 3.

#   There is also a site-specific covariate Pasture indicating whether the surrounding^
#   matrix is either predominately the modified habitat (farm pasture, Pasture =1) or
#   "native" grassland (tussock, Pasture = 0).

# 2018-08-18 Code submitted by Carl James Schwarz (cschwarz.stat.sfu.ca@gmail.com)

# Using the RPresence package

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

library(car)
library(readxl)
library(RPresence)
library(ggplot2)

# Get the RPResence additional functions 
source(file.path("..","..","..","AdditionalFunctions","Rpresence.additional.functions.R"))

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
sum(is.na(input.history))

site.covar <- data.frame(Site=1:nrow(input.data),
                         Pasture=car::recode(input.data$X__18,
                                            "1='Modified'; 0='Native';"))
head(site.covar)

# Five years with 3 visits. Don't need same number of visits/season
Nvisits.per.season  <- rep(3,5) # five years with 3 visits. Don't need same number of visits/season

# Create the *.pao file
skink.pao <- RPresence::createPao(input.history,
                                unitcov=site.covar,
                                nsurveyseason=Nvisits.per.season,
                                title='skink SSMS')
skink.pao

#-----
# Fit some models.
# We start with the first formulation in terms of psi, gamma, epsilon, p (do.1_)
# Note that formula DO NOT HAVE AN = SIGN
mod.fit <- RPresence::occMod(
           model=list(psi~1, gamma~1, epsilon~1, p~1), 
                        type="do.1", data=skink.pao)
summary(mod.fit)

mod.fit <- RPresence.add.derived(mod.fit)


# Estimate of initial occupance
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

