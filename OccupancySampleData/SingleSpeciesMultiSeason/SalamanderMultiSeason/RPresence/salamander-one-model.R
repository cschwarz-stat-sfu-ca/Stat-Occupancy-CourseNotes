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
library(RPresence)
library(ggplot2)

# Get the RPResence additional functions 
source(file.path("..","..","..","AdditionalFunctions","Rpresence.additional.functions.R"))

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





# Number of visits in each season
Nvisits.per.season  <- c(5,5) 

# Create the *.pao file
salamander.pao <- RPresence::createPao(input.history,
                                nsurveyseason=Nvisits.per.season,
                                unitcov=site.covar,
                                title='salamander SSMS')
salamander.pao

#-----
# Fit a models.
# We start with the first formulation in terms of psi, gamma, epsilon, p (do.1_)
# Note that formula DO NOT HAVE AN = SIGN
mod.fit<- RPresence::occMod(
           model=list(psi~1, gamma~1, epsilon~1, p~1), 
                                        type="do.1", data=salamander.pao)
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



