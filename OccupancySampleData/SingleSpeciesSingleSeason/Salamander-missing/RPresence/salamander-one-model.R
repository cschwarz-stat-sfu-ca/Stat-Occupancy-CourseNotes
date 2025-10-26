# Blue Ridge Salamander - with missing data
# Notice how the data are read in with the missing values.

# Single Species Single Season Occupancy

# Fitting a single model

# 2018-08-15 Code contributed by Carl James Schwarz (cschwarz.stat.sfu.cs@gmail.com)

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

library(readxl)
library(RPresence)
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


# Extract the history records
input.history <- input.data[, 1:5] # the history extracted
head(input.history)

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.history)
ncol(input.history)
range(input.history, na.rm=TRUE) # check that all values are either 0 or 1
sum(is.na(input.history))    # are there any missing values?


# Create the *.pao file
salamander.pao <- RPresence::createPao(input.history,
                                       title='Salamander SSSS')
salamander.pao



# Fit a model
# Note that formula DO NOT HAVE AN = SIGN
mod.fit <- RPresence::occMod(model=list(psi~1, p~1),
                              type="so", 
                              data=salamander.pao)
summary(mod.fit)

# Look the objects returned in more details
names(mod.fit)


# look at estimated occupancy probability. RPresence gives for EACH site in case it depends on covariates
mod.fit$beta$psi  # on the logit scale

mod.fit$real$psi  # on the regular 0-1 scale for each site
mod.fit$real$psi[1:5,]


# look at the estimated probability of detection. It gives an estimate for every site at very visit
mod.fit$real$p[1:5,]

# extract the detection probabilities for each visit. The row names have the unit names
mod.fit.p   <- mod.fit$real$p[grepl("unit1$", row.names(mod.fit$real$p)),]
mod.fit.p

# Look at the posterior probability of detection
names(mod.fit$derived)

mod.fit$derived$psi_c

# alternatively
RPresence::print_one_site_estimates(mod.fit, site = 1)

fitted(mod.fit, param="psi")[1:5,]
fitted(mod.fit, param="p")[1:5,]

