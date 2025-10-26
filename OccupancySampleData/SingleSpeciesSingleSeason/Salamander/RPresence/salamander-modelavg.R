# Single Species Single Season Occupancy 

# Salamander Example with model averaging with individual model fits

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
                                 sheet="CompleteData",
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



#-----
# Fit some models.
# Note that formula DO NOT HAVE AN = SIGN
mod.pdot <- RPresence::occMod(model=list(psi~1, p~1),
                              type="so", data=salamander.pao)
summary(mod.pdot)

names(mod.pdot)
mod.pdot$beta$psi

# look at estimated occupancy probability. RPresence gives for EACH site in case it depends on covariates
mod.pdot$real$psi[1:5,]

mod.pdot.psi <-mod.pdot$real$psi[1,]  # occupancy probability
mod.pdot.psi


# look at the estimated probability of detection. It gives an estimate for every site at very visit
head(mod.pdot$real$p)

mod.pdot.p   <- mod.pdot$real$p[grepl("unit1$", row.names(mod.pdot$real$p)),]
mod.pdot.p

# Look at the posterior probability of detection
names(mod.pdot$derived)

mod.pdot$derived$psi_c

# alternatively
RPresence::print_one_site_estimates(mod.pdot, site = 1)

#-------
# Model where p(t) varies across survey occasions
# 
mod.pt <- RPresence::occMod(model=list(psi~1, p~SURVEY), 
                            type="so", data=salamander.pao)
summary(mod.pt)

mod.pt$real$psi[1,]
mod.pt.p   <- mod.pt$real$p[grepl("unit1$", row.names(mod.pt$real$p)),]
mod.pt.p

print_one_site_estimates(mod.pt, site = 1)
fitted(mod.pt, param="psi")[1:5,]

#-----
# Fit a model with the detection probability equal in first 2 occasions and last 3 occasions.
# We need to define a survey covariate that has two levels 
# This covariate needs to be "stacked" so that sites1...site39 for survey occastion 1
# are then followed by covariate at survey occastion 2 for sites1...site39, etc

survey.cov <- data.frame(site=rep(1:nrow(input.history), ncol(input.history)),
                         visit=rep(1:ncol(input.history), each=nrow(input.history)),
                         d=rep( c("d1","d1","d2","d2","d2"),  each=nrow(input.history)))
head(survey.cov)

survey.cov[c(1:4, 37:41),]


mod.pcustom <- RPresence::occMod(model=list(psi~1, p~d), 
                      cov.list=list(p.cov=survey.cov), 
                      type="so", 
                      data=salamander.pao)

mod.pcustom$real$psi[1,]
mod.pcustom$real$p[grepl("unit1$", row.names(mod.pt$real$p)),]


#------
# Model averaging
models<-list(mod.pdot,mod.pt,mod.pcustom)
results<-RPresence::createAicTable(models)
summary(results)

RPresence::modAvg(results, param="psi")[1,]




