# Single Species Single Season Occupancy models using RPresence 

# Blue Gross Beaks.
#Downloaded from https://sites.google.com/site/asrworkshop/home/schedule/r-occupancy-1

#An occupancy study was made on Blue Grosbeaks (Guiraca caerulea) 
# on 41 old fields planted to longleaf pines (Pinus palustris) 
# in southern Georgia, USA. 

# Surveys were 500 m transects across each field 
# and were completed three times during the breeding season in 2001.

# Columns in the file are:
#    field - field number
#    v1, v2, v3 -  detection histories for each site on each of 3 visit during the 2001 breeding season.    
#    field.size - size of the files
#    bqi - Enrollment in bobwihte quail initiative; does occupancy increase if field belongs to this initiative?
#    crop.hist - crop history
#    crop1, crop2 - indicator variables for the crop history
#    count1, count2, count3 - are actual counts of birds detected in each visit

# 2018-08-15 Code contributed by Carl James Schwarz (cschwarz.stat.sfu.cs@gmail.com)



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

input.data <- read.csv(file.path("..","blgr.csv"), 
                       header=TRUE, as.is=TRUE, strip.white=TRUE) 
head(input.data)

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.data)
range(input.data[, c("v1","v2","v3")], na.rm=TRUE)
sum(is.na(input.data[, c("v1","v2","v3")]))

input.history <- input.data[, c("v1","v2","v3")]
head(input.history)

site.covar <- input.data[, c("field","field.size")]
site.covar$logFS <- log(site.covar$field.size)
head(site.covar)


# Create the *.pao file
grossbeak.pao <- RPresence::createPao(input.history,
                                      unitcov=site.covar,
                                      title='Grossbeak SSSS')
grossbeak.pao


#-----
# Fit some models.
# Note that formula DO NOT HAVE AN = SIGN
mod.pdot <- RPresence::occMod(model=list(psi~1, p~1),
                              type="so", data=grossbeak.pao)
summary(mod.pdot)

names(mod.pdot)
mod.pdot$beta$psi


# look at estimated occupancy probability. RPresence gives for EACH site in case it depends on covariates
head(mod.pdot$real$psi)

mod.pdot.psi <-mod.pdot$real$psi[1,]  # occupancy probability
mod.pdot.psi


# look at the estimated probability of detection. It gives an estimate for every site at very visit
head(mod.pdot$real$p)

mod.pdot.p   <- mod.pdot$real$p[seq(1, by=nrow(input.history), length.out=ncol(input.history)),]
mod.pdot.p

# Look at the posterior probability of detection
names(mod.pdot$derived)

mod.pdot$derived$psi_c

# alternatively
RPresence::print_one_site_estimates(mod.pdot, site = 1)

#-----
# Fit some models.
# Note that formula DO NOT HAVE AN = SIGN
mod.psilogFS.pdot <- RPresence::occMod(model=list(psi~logFS, p~1),
                              type="so", data=grossbeak.pao)
summary(mod.psilogFS.pdot)

names(mod.psilogFS.pdot)
mod.psilogFS.pdot$beta$psi


# look at estimated occupancy probability. RPresence gives for EACH site in case it depends on covariates
head(mod.psilogFS.pdot$real$psi)

mod.psilogFS.pdot.psi <-mod.psilogFS.pdot$real$psi[1,]  # occupancy probability
mod.psilogFS.pdot.psi

# plot of psi vs logfs #1
# individual psi
mod.psilogFS.pdot$real$psi
# covariate values
site.covar
both <- cbind(mod.psilogFS.pdot$real$psi, site.covar)

ggplot(data=both, aes(x=logFS, y=est))+
   geom_point()+
   geom_ribbon(aes(ymin=lower_0.95 , ymax=upper_0.95  ), alpha=0.2)


# plot #2
# what are beta values
mod.psilogFS.pdot$beta$psi
plotdata <- data.frame(logFS=seq(1,5,.1))
plotdata$logitpsi <-mod.psilogFS.pdot$beta$psi[1,1]+
         mod.psilogFS.pdot$beta$psi[2,1]*plotdata$logFS
plotdata$psi <- 1/(1+exp(-plotdata$logitpsi))
ggplot(data=plotdata, aes(x=logFS, y=psi))+
  geom_point()



# look at the estimated probability of detection. It gives an estimate for every site at very visit
head(mod.psilogFS.pdot$real$p)

mod.psilogFS.pdot.p   <- mod.psilogFS.pdot$real$p[seq(1, by=nrow(input.history), length.out=ncol(input.history)),]
mod.psilogFS.pdot.p

# Look at the posterior probability of detection
names(mod.psilogFS.pdot$derived)

mod.psilogFS.pdot$derived$psi_c

# alternatively
RPresence::print_one_site_estimates(mod.psilogFS.pdot, site = 1)

#-------
# Model where p(t) varies across survey occasions
# 
mod.pt <- RPresence::occMod(model=list(psi~1, p~SURVEY), type="so", data=grossbeak.pao)
summary(mod.pt)

mod.pt$real$psi[1,]
mod.pt$real$p[seq(1, by=nrow(input.history), length.out=ncol(input.history)),]

print_one_site_estimates(mod.pt, site = 1)
fitted(mod.pt, param="psi")

RPresence::print_one_site_estimates(mod.pt, site = 1)

#------
# Model averaging
models<-list(mod.pdot,mod.pt)
results<-RPresence::createAicTable(models)
summary(results)

RPresence::modAvg(results, param="psi")[1,]


