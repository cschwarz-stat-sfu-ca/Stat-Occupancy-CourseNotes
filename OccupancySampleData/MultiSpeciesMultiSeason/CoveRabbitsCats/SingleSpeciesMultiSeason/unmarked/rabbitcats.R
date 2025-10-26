# Single Species, Multi Season Occupancy analyais

#   Cove MV, Gardner B, Simons TR, O'Connell AF (2018) 
#   Co-occurrence dynamics of endangered Lower Keys marsh rabbits and free-ranging domestic cats: 
#   Prey responses to an exotic predator removal program. 
#   Ecology and Evolution 8(8): 4042-4052. 
#   https://doi.org/10.1002/ece3.3954

# Data obtained from the Dryad data package:
#
#  Cove MV, Gardner B, Simons TR, O'Connell AF (2018) 
#  Data from: Co-occurrence dynamics of endangered Lower Keys marsh rabbits and free-ranging domestic cats: 
#  prey responses to an exotic predator removal program. 
#  Dryad Digital Repository. 
#   https://doi.org/10.5061/dryad.748pd64

# Briefly: between 2013 and 2015, camera traps were used to survey marsh
# rabbits and free-ranging
# cats at 84 sites in the National Key Deer Refuge, Big Pine
# Key, Florida, USA. We used dynamic occupancy models to determine factors associated
# with marsh rabbit occurrence, colonization, extinction, and the co-occurrence
# of marsh rabbits and cats during a period of predator removal

# 2019-12-23 Code contributed by Carl James Schwarz (cschwarz.stat.sfu.cs@gmail.com)

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

#  load libraries 

library(readxl)
library(unmarked)
library(ggplot2)


# get the data read in
# Data for detections should be a data frame with rows corresponding to sites
# and columns to visits.
# The usual 1=detected; 0=not detected; NA=not visited is used.

input.data <- read.csv(file.path("..","..","Cove_etal_Ecology_Evolution_2018_data.csv"), 
                       header=TRUE, as.is=TRUE, strip.white=TRUE, na.string=".")
head(input.data)


#####################################################################
#####################################################################
#####################################################################
# First model dynamic occupancy for the rabbits
# See Table 1 of the paper

rabbit.detect.vars <- c("Rabbits.2013.2014","X", paste0("X.",1:14),"Rabbits.2015",paste0("X.",15:29))
rabbit.input.history <- input.data[, rabbit.detect.vars]


# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(rabbit.input.history)
ncol(rabbit.input.history)
range(rabbit.input.history, na.rm=TRUE)
sum(is.na(rabbit.input.history))
sum(rabbit.input.history, na.rm=TRUE)

# Site covariates.
# Groan, these must be 10 character or less (a limitation in RMark)
# we need to rename some of the variables
input.data <- plyr::rename(input.data, c("salt_button"="salt_butt",
                                         "CatIndividuals_2013_14"="icats1",
                                          "dist_cat2014_z"="dcat2014z",
                                           "bin_2014cat"="icat2014"))
names(input.data)

unique(input.data$DESCRIPTION)
xtabs(~DESCRIPTION+salt_butt, data=input.data,exclude=NULL, na.action=na.pass)

global.covar.vars <- c("HumanTrail","Rabbitat","dist_dev_z","freshwater","salt_butt","hectares_z","LiDar_z")
setdiff(global.covar.vars, names(input.data))

input.data$CatIndividuals_2013_14
input.data$CatIndividuals_2015
input.data$delta.icats <- input.data$CatIndividuals_2013_14 - input.data$CatIndividuals_2015
input.data$delta.icats


input.data$Total.CatCaps_2013_14
input.data$Total.CatCaps_2015

# distance to sites with cat removed in 2014 vs binary indicator
input.data$icat2014
input.data$dcat2014z
ggplot(data=input.data, aes(x=dcat2014z, y=icat2014))+
   ggtitle("indicatior if cats captures within 500m buffer of camera trap")+
   geom_point()

other.covar.vars <- c("icats1",    # number of individual cats in year 1
                      "dcat2014z", # distance to site with a cat in 2014 
                      "icat2014")  # indicator if cat captured within 500m of the site in 2014
setdiff(other.covar.vars, names(input.data))

site.covar <- input.data[, c(global.covar.vars,other.covar.vars)]
row.names(site.covar) <- NULL
head(site.covar)


# Create the "season" matrix so that you know which history is in which season
year.covar <- matrix(c("1","2"),nrow=nrow(rabbit.input.history), ncol=2, byrow=TRUE)
head(year.covar)

rabbit.data <-unmarkedMultFrame(y=rabbit.input.history,
                                siteCovs=site.covar, 
                                yearlySiteCovs=as.data.frame(year.covar),
                                numPrimary=2)   # two primary sessions each with 16 secondary sessions


## Fit some models using unmarked

fm1 <- colext(psiformula=~1,
              gammaformula=~1,
              epsilonformula=~1,
              pformula=~1,
              data=rabbit.data)
fm1

fm2 <- colext(psiformula=~~HumanTrail+Rabbitat+dist_dev_z+freshwater+salt_butt+hectares_z+LiDar_z+icats1,
              gammaformula=~1,
              epsilonformula=~1,
              pformula=~HumanTrail,
              data=rabbit.data)
fm2


#back transformations
backTransform(fm1,'det')
backTransform(fm1,"psi")
backTransform(fm1,"col")
backTransform(fm1,"ext")

