# Single Species Single Season Occupancy 

# Weta Example

#  RPresence package

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

input.history <- readxl::read_excel(file.path("..","weta.xls"),
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
site_covar <- readxl::read_excel(file.path("..","weta.xls"),
                                 sheet="site_covar",
                                 na="-",
                                 col_names=TRUE)  # notice col_names in row 1 of table. 
head(site_covar)

# Create an alternate site level covariate that is a categorical variable rather 
# than indicator variables
site_covar$BrowCat <- site_covar$BrowseCat
xtabs(~BrowCat, data=site_covar,exclude=NULL, na.action=na.pass)

head(site_covar)


# Get the individual covariates. 
obs1 <- readxl::read_excel(file.path("..","weta.xls"),
                           sheet="Obs1",
                           na="-",
                           col_names=FALSE) 
obs2 <- readxl::read_excel(file.path("..","weta.xls"),
                           sheet="Obs2",
                           na="-",
                           col_names=FALSE) 
obs3 <- readxl::read_excel(file.path("..","weta.xls"),
                           sheet="Obs3",
                           na="-",
                           col_names=FALSE) 

Obs <- obs1*1 + obs2*2 + obs3*3
head(Obs)

# Observational covariate needs to be "stacked" so that sites1...siteS for survey occastion 1
# are then followed by covariate at survey occastion 2 for sites1...siteS, etc

survey.cov <- data.frame(site=rep(1:nrow(input.history) , ncol(input.history)),
                         visit=as.character(rep(1:ncol(input.history), each=nrow(input.history))),  # notice we make a character 
                         obs1 =as.vector(unlist(obs1)),
                         obs2 =as.vector(unlist(obs2)),
                         obs3 =as.vector(unlist(obs3)),
                         Obs  =paste("O",as.vector(unlist(Obs)),sep=""),    # notice we make a character string
                         stringsAsFactors=FALSE)
survey.cov$Obs[ grepl("NA",survey.cov$Obs)] <- NA
head(survey.cov[,c("visit","site","Obs")],n=100)
str(survey.cov)

# check that missing values in history and observer covariates align
select <- is.na(as.vector(unlist(input.history)))
survey.cov[select,]
sum(is.na(survey.cov[!select,]))


# Create the *.pao file
weta.pao <- RPresence::createPao(input.history,
                                 unitcov=site_covar,
                                 survcov=survey.cov,
                                 title='weta SSSS')
weta.pao

#--------------------------------------------------------------------------
# Fit some models.
# Note that formula DO NOT HAVE AN = SIGN
mod.pdot <- RPresence::occMod(model=list(psi~1, p~1), type="so", data=weta.pao)
summary(mod.pdot)

# look at estimated occupancy probability. RPresence gives for EACH site in case it depends on covariates
mod.pdot.psi <-mod.pdot$real$psi[1,]  # occupancy probability
mod.pdot.psi


# look at the estimated probability of detection. It gives an estimate for every site at very visit
mod.pdot.p   <- mod.pdot$real$p[seq(1, by=nrow(input.history), length.out=ncol(input.history)),]
mod.pdot.p

# alternatively
RPresence::print_one_site_estimates(mod.pdot, site = 1)


#--------------------------------------------------------------------------
# Impact of browse.
# Use the two indicator variables or the categorized column
# This give identical results

# Cell means models.
# We need to EXCLUDE the intercept if we are using the two indicator variables.
mod.pdot.psiB.1 <- RPresence::occMod(model=list(psi~-1+Browsed+Unbrowsed, p~1), type="so", data=weta.pao)
summary(mod.pdot.psiB.1)

names(mod.pdot.psiB.1)

mod.pdot.psiB.1.psi <-mod.pdot.psiB.1$real$psi  # occupancy probability
mod.pdot.psiB.1.psi[1:5,]

# This is the logit occupancy for each browse category
mod.pdot.psiB.1$dmat$psi

mod.pdot.psiB.1$beta$psi
mod.pdot.psiB.1.oddsratio.browse <- exp( sum(c(1,-1)*mod.pdot.psiB.1$beta$psi))
mod.pdot.psiB.1.oddsratio.browse

# Alternate way for cell means models.
mod.pdot.psiB.4 <- RPresence::occMod(model=list(psi~-1+BrowCat, p~1), type="so", data=weta.pao)
summary(mod.pdot.psiB.4)

names(mod.pdot.psiB.4)

mod.pdot.psiB.4.psi <-mod.pdot.psiB.4$real$psi  # occupancy probability
mod.pdot.psiB.4.psi[1:5,]

# This is the logit occupancy for each browse category
mod.pdot.psiB.4$dmat$psi

mod.pdot.psiB.4$beta$psi
mod.pdot.psiB.4.oddsratio.browse <- exp( sum(c(1,-1)*mod.pdot.psiB.4$beta$psi))
mod.pdot.psiB.4.oddsratio.browse



# We can use just a single Browser indicator which is the cell effects model.
mod.pdot.psiB.2 <- RPresence::occMod(model=list(psi~Browsed, p~1), type="so", data=weta.pao)
summary(mod.pdot.psiB.2)

names(mod.pdot.psiB.2)

mod.pdot.psiB.2.psi <-mod.pdot.psiB.2$real$psi  # occupancy probability
mod.pdot.psiB.2.psi[1:5,]

# This is the logit occupancy the baseline and the log(odds ratio) of the two 
# browse classes
mod.pdot.psiB.2$dmat$psi

mod.pdot.psiB.2$beta$psi
mod.pdot.psiB.2.oddsratio.browse <- exp( sum(c(0,1)*mod.pdot.psiB.2$beta$psi))
mod.pdot.psiB.2.oddsratio.browse



# We can use just a categorical browser variable (preferred).
mod.pdot.psiB.3 <- RPresence::occMod(model=list(psi~BrowCat, p~1), type="so", data=weta.pao)
summary(mod.pdot.psiB.3)

names(mod.pdot.psiB.3)

mod.pdot.psiB.3.psi <-mod.pdot.psiB.3$real$psi  # occupancy probability
mod.pdot.psiB.3.psi[1:5,]

# This is the logit occupancy the baseline and the log(odds ratio) of the two 
# browse classes
mod.pdot.psiB.3$dmat$psi
mod.pdot.psiB.3$beta$psi
mod.pdot.psiB.3.oddsratio.browse <- exp( sum(c(0,-1)*mod.pdot.psiB.3$beta$psi))
mod.pdot.psiB.3.oddsratio.browse


#----------------------
# Impact of browse on detectability as well?.

mod.pB.psiB <- RPresence::occMod(model=list(psi~BrowCat, p~BrowCat), type="so", data=weta.pao)
summary(mod.pB.psiB)

names(mod.pB.psiB)

mod.pB.psiB.psi <-mod.pB.psiB$real$psi  # occupancy probability
mod.pB.psiB.psi[1:5,]

mod.pB.psiB.p <-mod.pB.psiB$real$p  # detection probability
mod.pB.psiB.p[1:5,]


# This is the logit occupancy for each browse category
mod.pB.psiB$dmat$psi
mod.pB.psiB$beta$psi
mod.pB.psiB.oddsratio.browse <- exp( sum(c(0,-1)*mod.pB.psiB$beta$psi))
mod.pB.psiB.oddsratio.browse


#-------
# Model where p(t) varies across survey occasions but psi(browse) 
# Not sure why SURVEY key work no longer works for time varying factor
# but we have our own covariate (visit). Be sure to declare "visit" as a character or as a actor
# 
mod.pt.psiB <- RPresence::occMod(model=list(psi~BrowCat, p~visit), type="so", data=weta.pao)
summary(mod.pt.psiB)

mod.pt.psiB$real$psi[1:5,]
mod.pt.psiB$real$p[seq(1, by=nrow(input.history), length.out=ncol(input.history)),]

print_one_site_estimates(mod.pt.psiB, site = 1)
fitted(mod.pt.psiB, param="psi")[1,]



#-------
# Model where p varies by observer but constant over time psi(browse) 
# 
mod.pO.psiB <- RPresence::occMod(model=list(psi~BrowCat, p~Obs), type="so", data=weta.pao)
summary(mod.pt.psiB)

mod.pO.psiB$real$psi[1:5,]
mod.pO.psiB$real$p[seq(1, by=nrow(input.history), length.out=ncol(input.history)),]

print_one_site_estimates(mod.pO.psiB, site = 1)
fitted(mod.pO.psiB, param="psi")[1:5,]


#-------
# Model where p varies by observer + visit but constant over time psi(browse)
# besure that obs and visit are declared as character
# 
mod.pOpV.psiB <- RPresence::occMod(model=list(psi~BrowCat, p~Obs+visit), type="so", data=weta.pao)
summary(mod.pOpV.psiB)

mod.pOpV.psiB$real$psi[1:5,]
mod.pOpV.psiB$real$p[seq(1, by=nrow(input.history), length.out=ncol(input.history)),]

print_one_site_estimates(mod.pOpV.psiB, site = 1)
fitted(mod.pOpV.psiB, param="psi")[1:5,]

#-------
# Other models 

mod.pOpVpB.psiB <- RPresence::occMod(model=list(psi~BrowCat, p~Obs+visit+BrowCat), type="so", data=weta.pao)

mod.pOpVpB.psi. <- RPresence::occMod(model=list(psi~1, p~Obs+visit+BrowCat), type="so", data=weta.pao)

mod.pOpVpB.psi. <- RPresence::occMod(model=list(psi~1, p~Obs+visit), type="so", data=weta.pao)



#------
# Model averaging
models<-list(mod.pdot, 
             mod.pdot.psiB.1,
             mod.pB.psiB, 
             mod.pt.psiB,
             mod.pO.psiB,
             mod.pOpV.psiB,
             mod.pOpVpB.psiB,
             mod.pOpVpB.psi.,
             mod.pOpVpB.psi.
)
results<-RPresence::createAicTable(models)
summary(results)

RPresence::modAvg(results, param="psi")[1:5,]

