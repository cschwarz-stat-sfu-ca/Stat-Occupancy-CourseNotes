# American Toad

# Extracted from
#    Darryl I. MacKenzie, et al. 2002.
#    Estimating site occupancy rates when detection probabilities are less than one.
#    Ecology 83:2248-2255.
#    doi:10.1890/0012-9658(2002)083[2248:ESORWD]2.0.CO;2]

# 29 sites with 82 sampling occasions in 2000.
# Volunteers visited sites and recorded presence/absence of toads by calls.
# Habitat (type of pond, permanent or ephemeral) and temperature at visit recorded.

# Single Species Single Season Occupancy

# Fitting a single model using RPresence

# 2018-08-15 Code contributed by Carl James Schwarz (cschwarz.stat.sfu.cs@gmail.com)

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

library(car)
library(readxl)
library(RPresence)
library(ggplot2)

# get the data read in
# Data for detections should be a data frame with rows corresponding to sites
# and columns to visits.
# The usual 1=detected; 0=not detected; NA=not visited is used.

input.data <- readxl::read_excel("../AmericanToad.xls",
                                 sheet="AmToadDetectionHistories",
                                 na="-",
                                 col_names=FALSE)  # notice no column names in row 1 of data file. 

head(input.data)


# Extract the history records
input.history <- input.data # the history extracted
head(input.history)

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.history)
ncol(input.history)
range(input.history, na.rm=TRUE) # check that all values are either 0 or 1
sum(is.na(input.history))    # are there any missing values?


# Get the pond information
pond.data <- readxl::read_excel("../AmericanToad.xls",
                                sheet="Pond",
                                na="-",
                                col_names=TRUE)  
head(pond.data)
pond.data$Pond <- as.factor(car::recode(pond.data$Pond,
                              "1='P'; 0='E'; " ))
head(pond.data)


# Get the temperature data
temp.data <- readxl::read_excel("../AmericanToad.xls",
                                sheet="AmToadTemperature",
                                na="-",
                                col_names=FALSE)  # notice no column names in row 1 of data file. 
head(temp.data)
temp.data[ is.na(input.history)] <- NA

# Convert to a survey covariate. You need to stack the data by columns
survey.cov <- data.frame(Site=1:nrow(temp.data),
                         visit=as.factor(rep(1:ncol(temp.data), each=nrow(temp.data))),
                         Temperature =unlist(temp.data), stringsAsFactors=FALSE)

head(survey.cov)


# Create the *.pao file
amtoad.pao <- RPresence::createPao(input.history,
                                   unitcov=pond.data,
                                   survcov=survey.cov,
                                  title='American Toad SSSS')
amtoad.pao



# Fit a model
# Note that formula DO NOT HAVE AN = SIGN
mod.fit <- RPresence::occMod(model=list(psi~Pond, p~Temperature),
                              type="so", 
                              data=amtoad.pao)
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

