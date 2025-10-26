# Single Species Single Season Occupancy models using RPresence and unmarked packages and JAGS

# Weta Example

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

#  RPresence package

library(readxl)
library(RPresence)
library(ggplot2)

# get the data read in
# Data for detections should be a data frame with rows corresponding to sites
# and columns to visits.
# The usual 1=detected; 0=not detected; NA=not visited is used.

input.history <- readxl::read_excel("Weta_pg116.xls",
                                 sheet="detection_histories",
                                 na="-",
                                 col_names=FALSE)  # notice no column names in row 1 of data file. 

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.history)
ncol(input.history)
range(input.history, na.rm=TRUE)
sum(is.na(input.history))

head(input.history)

# Get the site level covariates
site_covar <- readxl::read_excel("Weta_pg116.xls",
                                 sheet="site_covar",
                                 na="-",
                                 col_names=TRUE)  # notice col_names in row 1 of table. 


# Create an alternate site level covariate that is a categorical variable rather 
# than indicator variables
site_covar$BrowCat <- paste(c("","B")[1+unlist(site_covar[,1])], c("","N")[1+unlist(site_covar[,2])], sep="")
xtabs(~BrowCat, data=site_covar,exclude=NULL, na.action=na.pass)
colSums(site_covar[,1:2])

head(site_covar)


# Get the individual covariates. 
obs1 <- readxl::read_excel("Weta_pg116.xls",
                                 sheet="Obs1",
                                 na="-",
                                 col_names=FALSE) 
obs2 <- readxl::read_excel("Weta_pg116.xls",
                                 sheet="Obs2",
                                 na="-",
                                 col_names=FALSE) 
obs3 <- readxl::read_excel("Weta_pg116.xls",
                                 sheet="Obs3",
                                 na="-",
                                 col_names=FALSE) 

obs <- obs1*1 + obs2*2 + obs3*3
head(obs)

# Observational covariate needs to be "stacked" so that sites1...siteS for survey occastion 1
# are then followed by covariate at survey occastion 2 for sites1...siteS, etc

survey.cov <- data.frame(site=rep(1:nrow(input.history) , ncol(input.history)),
                         visit=rep(1:ncol(input.history), each=nrow(input.history)),
                         obs1 =as.vector(unlist(obs1)),
                         obs2 =as.vector(unlist(obs2)),
                         obs3 =as.vector(unlist(obs3)),
                         obs  =as.character(as.vector(unlist(obs))), stringsAsFactors=FALSE)
head(survey.cov)

# check that missing values in history and observer covariates align
select <- is.na(as.vector(unlist(input.history)))
survey.cov[select,]
sum(is.na(survey.cov[!select,]))

weta.pao.nas <-RPresence::createPao(input.history,
                                       unitcov=site_covar,
                                       survcov=survey.cov,
                                       title='weta SSSS')
weta.pao.nas

# The missing values in the survey covariates must be filled with dummy
# values to avoid problems in fitting the models that depend on them
survey.cov[ is.na(survey.cov)] <- -1
survey.cov[select,]
sum(is.na(survey.cov[!select,]==-1))


# Create the *.pao file
weta.pao <- RPresence::createPao(input.history,
                                       unitcov=site_covar,
                                       survcov=survey.cov,
                                       title='weta SSSS')
weta.pao

weta.pao.nocov <- RPresence::createPao(input.history,
                                       title='weta SSSS')
weta.pao.nocov


##### Problem #1 = cannot use SURVEY in time dependent model if you have covariates

mod.pdot.withcov <- RPresence::occMod(model=list(psi~1, p~SURVEY), type="so", data=weta.pao)

mod.pdot.nocov   <- RPresence::occMod(model=list(psi~1, p~SURVEY), type="so", data=weta.pao.nocov)



##### Probem #2 - missing values in covariates cause problems 

mod.pO.psiB <- RPresence::occMod(model=list(psi~BrowCat, p~obs), type="so", data=weta.pao.nas)

mod.pO.psiB <- RPresence::occMod(model=list(psi~BrowCat, p~obs), type="so", data=weta.pao)


