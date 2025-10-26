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
library(RPresence)
library(ggplot2)

# Get the Mark additional functions 
source(file.path("..","..","..","..","AdditionalFunctions","RPresence.additional.functions.R"))

# get the data read in
# Data for detections should be a data frame with rows corresponding to sites
# and columns to visits.
# The usual 1=detected; 0=not detected; NA=not visited is used.

input.data <- read.csv(file.path("..","..","Cove_etal_Ecology_Evolution_2018_data.csv"), 
                       header=TRUE, as.is=TRUE, strip.white=TRUE, na.string=".")
head(input.data)

.
# Groan, covariate names must be 10 characters or less (a limitation in RMark)
# we need to rename some of the variables
input.data <- plyr::rename(input.data, c("salt_button"           ="salt_butt",
                                         "CatIndividuals_2013_14"="icats1",
                                          "dist_cat2014_z"       ="dcat2014z",
                                           "bin_2014cat"         ="icat2014"))
names(input.data)

# the paper collapses the different habitats into 3 categories, and creates indicator variables
unique(input.data$DESCRIPTION)
xtabs(~DESCRIPTION+salt_butt,  data=input.data,exclude=NULL, na.action=na.pass)
xtabs(~DESCRIPTION+freshwater, data=input.data,exclude=NULL, na.action=na.pass)
xtabs(~DESCRIPTION+upland,     data=input.data,exclude=NULL, na.action=na.pass)



# compute the difference in individual cats. These are problematic as the values are the same
# in the two sessions
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


# The author defines a set of global covariates that are used
global.covar.vars <- c("HumanTrail","Rabbitat","dist_dev_z","freshwater","salt_butt","hectares_z","LiDar_z")
setdiff(global.covar.vars, names(input.data))

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


other.covar.vars <- c("icats1",    # number of individual cats in year 1
                      "dcat2014z", # distance to site with a cat in 2014 
                      "icat2014")  # indicator if cat captured within 500m of the site in 2014
setdiff(other.covar.vars, names(input.data))

site.covar <- input.data[, c(global.covar.vars,other.covar.vars)]
row.names(site.covar) <- NULL
head(site.covar)


#Create the RPresence data structure
# 2 years with 15 camera trapping days/year
max.visit.per.year <- 16
n.season <- 2

rabbit.pao <- RPresence::createPao(data=rabbit.input.history, 
                                unitcov=site.covar,
                                nsurveyseason=rep(max.visit.per.year,n.season))
rabbit.pao


# The paper claims there are 10 observed rabbit colonization events and 7 observed extinction events 
# Is this true? 
dim(rabbit.input.history)
detect.rabbit <- data.frame(detect.r1=apply(rabbit.input.history[, 1:16],1,max, na.rm=TRUE),
                            detect.r2=apply(rabbit.input.history[,17:32],1,max, na.rm=TRUE))
detect.rabbit$ch <- apply(detect.rabbit,1,paste0,collapse="")
head(detect.rabbit)
xtabs(~ch, data=detect.rabbit)


# all time-specific covariates must be added to the ddl in the loop below



# The top model is psi(global), gamma(ind change in cats) epsilon (dot) p(humtrail)
# If you look at the supplemental table for model 20 you see that
# this is psi(humtrail+rabbitat+devel+fresh+salt_butt+patch+LiDar)

# Notice the commas between the column and the placement of the quotes
# Define the models.
#    model.type can be RDOccupPE, RDOccupPG, do.1~time
#The random occupancy model cannot be fit here. See the other code in this directory.
# The reserved keywords can be found by looking at the ddl structures earlier

global <- paste0("~",paste0(global.covar.vars, collapse="+"),collapse="")

model.list.csv <- textConnection("
p,               psi,               gamma,      epsilon,  model.type
~HumanTrail,   ~global,                ~1,           ~1,       do.1
~HumanTrail,   ~global+icats1,         ~1,           ~1,       do.1
~HumanTrail,   ~global+icats1,         ~dcat2014z,   ~1,       do.1
~HumanTrail,   ~global+icats1,         ~icat2014,    ~1,       do.1
~1,            ~1,                     ~1,           ~1,       do.1")
  

model.list <- read.csv(model.list.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)
model.list$model.number < -1:nrow(model.list)
model.list

# substitute for some of the shortcuts
model.list$psi <- gsub("~global",global, model.list$psi, )
model.list


# fit the models
model.fits <- plyr::alply(model.list[2,], 1, function(x,detect.pao){
    cat("\n\n***** Starting ", unlist(x), "\n")
    if(x$model.type == 'do.1'){
       fit <- RPresence::occMod(model=list(as.formula(paste("psi",x$psi)),
                                        as.formula(paste("p"  ,x$p  )),
                                        as.formula(paste("gamma",x$gamma)),
                                        as.formula(paste("epsilon",x$epsilon))),
              data=detect.pao,type="do.1")
    }
    if(x$model.type == 'do.4'){
       fit <- RPresence::occMod(model=list(as.formula(paste("psi",x$psi)),
                                        as.formula(paste("p"  ,x$p  ))),
              data=detect.pao,type="do.4")
    }
    fit <- RPresence.add.derived(fit)
    fit

},detect.pao=rabbit.pao)


# Look at output from a specified model
model.number <- 2


names(model.fits[[model.number]])
names(model.fits[[model.number]]$real)
model.fits[[model.number]]$beta
names(model.fits[[model.number]]$derived)
model.fits[[model.number]]$derived$psi[1:10,]
model.fits[[model.number]]$real$gamma[1:5,]
model.fits[[model.number]]$real$epsilon[1:5,]

# Estimate of initial occupance
model.fits[[model.number]]$real$psi[grepl('unit1_', row.names(model.fits[[model.number]]$real$psi)),]

# Derived parameters - estimated occupancy for each unit in years 2....
names(model.fits[[model.number]]$derived)
model.fits[[model.number]]$derived$psi[ grepl('unit1_', row.names(model.fits[[model.number]]$derived$psi)),]

# Derived parameters - all of the psi stacked together
model.fits[[model.number]]$derived$all_psi[ grepl('unit1_', row.names(model.fits[[model.number]]$derived$all_psi)),]

# Estimate of  local extinction probability for each unit
model.fits[[model.number]]$real$epsilon[ seq(1, by=nrow(input.history), length.out=length(Nvisits.per.season)-1),]

# Estimate of  local colonization probability for each unit
model.fits[[model.number]]$real$gamma[ seq(1, by=nrow(input.history), length.out=length(Nvisits.per.season)-1),]

# Estimate of probability of detection at each time point for each unit
model.fits[[model.number]]$real$p[ grepl('unit1_', row.names(model.fits[[model.number]]$real$p), fixed=TRUE),]


# Get the change in occupancy
# Not yet possible to estimate the se of these values. May have to use bootstrapping.
model.fits[[model.number]]$derived$lambda [grepl('unit1_', row.names(model.fits[[model.number]]$derived$lambda),  fixed=TRUE),]
model.fits[[model.number]]$derived$lambdap[grepl('unit1_', row.names(model.fits[[model.number]]$derived$lambdap), fixed=TRUE),]



# collect models and make AIC table
aic.table <- RPresence::createAicTable(model.fits)
aic.table$table


# model averaging in the usual way

# initial occupancy
RPresence::modAvg(aic.table, param="psi")[1:5,]

# model averaging of derived parameters such as the occupancy at each time step
ma_all_psi <- RPresence.modAvg.derived(aic.table, param="all_psi")
ma_all_psi[grepl('unit1_', row.names(ma_all_psi),  fixed=TRUE),]

# Unable to get model averaged estimates of lambda because se is currently unknown


#####################################################################
#####################################################################
#####################################################################
# Next model dynamic occupancy for the cats
# See Tables S2 and S3 of the paper

cat.detect.vars <- c("Cats.2013.2014", paste0("X.",31:45),"Cats.2015",paste0("X.",46:60))
cat.input.history <- input.data[, cat.detect.vars]

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(cat.input.history)
ncol(cat.input.history)
range(cat.input.history, na.rm=TRUE)
sum(is.na(cat.input.history))
sum(cat.input.history, na.rm=TRUE)


other.covar.vars <- c("year_2013")  # indicator if site first measured in 2013 or 2014
setdiff(other.covar.vars, names(input.data))

site.covar <- input.data[, c(global.covar.vars,other.covar.vars)]
row.names(site.covar) <- NULL
head(site.covar)


cat.pao <- RPresence::createPao(data=cat.input.history, 
                                unitcov=site.covar,
                                nsurveyseason=rep(max.visit.per.year,n.season))
cat.pao


# The paper claims there are 9 observed cat colonization events and 27 observed extinction events 
# Is this true? 
dim(cat.input.history)
detect.cat <- data.frame(detect.r1=apply(cat.input.history[, 1:16],1,max, na.rm=TRUE),
                            detect.r2=apply(cat.input.history[,17:32],1,max, na.rm=TRUE))
detect.cat$ch <- apply(detect.cat,1,paste0,collapse="")
head(detect.cat)
xtabs(~ch, data=detect.cat)


# The top model is psi(global), gamma(ind change in cats) epsilon (dot) p(humtrail)
# If you look at the supplemental table for model 20 you see that
# this is psi(humtrail+catat+devel+fresh+salt_butt+patch+LiDar)

# Notice the commas between the column and the placement of the quotes
# Define the models.
#    model.type can be RDOccupPE, RDOccupPG, do.1~time
#The random occupancy model cannot be fit here. See the other code in this directory.
# The reserved keywords can be found by looking at the ddl structures earlier

global <- paste0("~",paste0(global.covar.vars, collapse="+"),collapse="")

model.list.csv <- textConnection("
p,                                 psi,               gamma,      epsilon,      model.type
~HumanTrail+Rabbitat+year_2013,   ~global,                ~1,      ~year_2013,       do.1")


model.list <- read.csv(model.list.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)
model.list$model.number <- 1:nrow(model.list)
model.list

# substitute for some of the shortcuts
model.list$psi <- gsub("~global",global, model.list$psi, )
model.list

# fit the models
model.fits <- plyr::alply(model.list, 1, function(x,detect.pao){
    cat("\n\n***** Starting ", unlist(x), "\n")
    if(x$model.type == 'do.1'){
       fit <- RPresence::occMod(model=list(as.formula(paste("psi",x$psi)),
                                        as.formula(paste("p"  ,x$p  )),
                                        as.formula(paste("gamma",x$gamma)),
                                        as.formula(paste("epsilon",x$epsilon))),
              data=detect.pao,type="do.1")
    }
    if(x$model.type == 'do.4'){
       fit <- RPresence::occMod(model=list(as.formula(paste("psi",x$psi)),
                                        as.formula(paste("p"  ,x$p  ))),
              data=detect.pao,type="do.4")
    }
    fit <- RPresence.add.derived(fit)
    fit

},detect.pao=cat.pao)


# Look at output from a specified model
model.number <- 1


names(model.fits[[model.number]])
names(model.fits[[model.number]]$real)
model.fits[[model.number]]$beta
names(model.fits[[model.number]]$derived)
model.fits[[model.number]]$derived$psi[1:10,]
model.fits[[model.number]]$real$gamma[1:5,]
model.fits[[model.number]]$real$epsilon[1:5,]

# Estimate of initial occupance
model.fits[[model.number]]$real$psi[grepl('unit1_', row.names(model.fits[[model.number]]$real$psi)),]

# Derived parameters - estimated occupancy for each unit in years 2....
names(model.fits[[model.number]]$derived)
model.fits[[model.number]]$derived$psi[ grepl('unit1_', row.names(model.fits[[model.number]]$derived$psi)),]

# Derived parameters - all of the psi stacked together
model.fits[[model.number]]$derived$all_psi[ grepl('unit1_', row.names(model.fits[[model.number]]$derived$all_psi)),]

# Estimate of  local extinction probability for each unit
model.fits[[model.number]]$real$epsilon[ seq(1, by=nrow(input.history), length.out=length(Nvisits.per.season)-1),]

# Estimate of  local colonization probability for each unit
model.fits[[model.number]]$real$gamma[ seq(1, by=nrow(input.history), length.out=length(Nvisits.per.season)-1),]

# Estimate of probability of detection at each time point for each unit
model.fits[[model.number]]$real$p[ grepl('unit1_', row.names(model.fits[[model.number]]$real$p), fixed=TRUE),]


# Get the change in occupancy
# Not yet possible to estimate the se of these values. May have to use bootstrapping.
model.fits[[model.number]]$derived$lambda [grepl('unit1_', row.names(model.fits[[model.number]]$derived$lambda),  fixed=TRUE),]
model.fits[[model.number]]$derived$lambdap[grepl('unit1_', row.names(model.fits[[model.number]]$derived$lambdap), fixed=TRUE),]



# collect models and make AIC table
aic.table <- RPresence::createAicTable(model.fits)
aic.table$table


# model averaging in the usual way

# initial occupancy
RPresence::modAvg(aic.table, param="psi")[1:5,]

# model averaging of derived parameters such as the occupancy at each time step
ma_all_psi <- RPresence.modAvg.derived(aic.table, param="all_psi")
ma_all_psi[grepl('unit1_', row.names(ma_all_psi),  fixed=TRUE),]

#####################################################################
#####################################################################
#####################################################################
# Finally  model dynamic occupancy for both species
# See Table 3 in the paper

# RPresence wants a 1 digit code

both.history <- 2*rabbit.input.history + 1*cat.input.history
colnames(both.history) <- NULL
both.history[1:2,]

input.chistory <- data.frame(lapply(as.data.frame(input.history2), as.character), stringsAsFactors=FALSE)


other.covar.vars <- c("year_2013")  # indicator if site first measured in 2013 or 2014
setdiff(other.covar.vars, names(input.data))

site.covar <- input.data[, c(global.covar.vars,other.covar.vars)]
row.names(site.covar) <- NULL
head(site.covar)

#Create the RMark data structure
# 2 years with 15 camera trapping days/year
max.visit.per.year <- 16
n.season <- 2

both.pao <- RPresence::createPao(data=both.history, 
                                unitcov=site.covar,
                                nsurveyseason=rep(max.visit.per.year,n.season))
both.pao



#-----
# Use the psiAB parameterization 
#    Parameters of the model are 
#        PsiA - Initial probability that site is occupied by A only
#        PsiBA- Initial probability that site is occupied by B if A is present
#        PsiBa- Initial probability that site is occupied by B is A is absent
#
#    If species are independent thatn psiBA = psiBa.
#
#        GammaA  - Colonization by A given B absent at t-1
#        GammaB  - Colonization by B given A absent at t-1
#        GammaAB - Colonization by A given B present at t-1
#        GammaBA - Colonization by B given A present at t-1
 
#        EpsilonA  - extinction of A given B absent at t-1
#        EpsilonB  - extinction of B given A absent at t-1
#        EpsilonAB - extinction of A given B present at t-1
#        EpsilonBA - extinction of B given A present at t-1
#
#    This give the transition matrx (see (2) in the paper but notice that the 
#    From/To sides of the matrix have been transposed below
#                           To State
#                 U                            A                    B                          AB
#    From   U  (1-GammaA)(1-GammaB)      GammaA(1-GammaB)    (1-GammaA)GammaB             GammaA GammaB
#    State  A   EpsilonA (1-GammaBA) (1-EpsilonA)(1-GammaBA)   EpsilonGammaBA         (1-EpsilonA)(GammaBA)
#           B  (1-gammaAB)EpsilonB   GammaAB (EpsilonB)      (1-GammaAB)(1-EpsilonB)   GammaAB(1-EpsilonB)
#          AB   EpsilonAB EpsilonBA (1-EpsilonAB)EpsilonBA    EpsilonAB(1-EpsionBA)     (1-EpsionAB)(1-EpsilonBA)
#
#
#    Detection parameters
#        pA    - probability of detection if A is alone in the site
#        pB    - probability of detection if B is alone in the site
#        rA    - probability of detecting A only given both on the site
#        rBA   - probability of detecting B given that A was detected and both are on the site
#        rBa   - probability of detecting B given that A not detected and both are on the site 
#    Ifspecies do not interact, then
#        rBA = rBa
#    and
#        rho = odds(rBA)/odds(rBa) = 1




# Get the list of models. NOtice NO equal signs here
# Here A=cat, B=rabbit
global <- paste0("~",paste0(global.covar.vars, collapse="+"),collapse="")

model.list.csv <- textConnection("
p,                                 psi,               gamma,      epsilon,      model.type
~1,   ~1,                ~1,      ~1,       do.2sp.1
~1,   ~SP+INT,                ~1,      ~1,       do.2sp.1
~1,   ~fw+dh+ps+SP+INT,                ~1,      ~1,       do.2sp.1")

model.list <- read.csv(model.list.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)
model.list$model.number <- 1:nrow(model.list)

model.list$psi  <- gsub("fw","freshwater",model.list$psi)
model.list$psi  <- gsub("ps","hectares_z",model.list$psi)  # patch size
model.list$psi  <- gsub("dh","dist_dev_z",model.list$psi)  # distance to development
model.list

# fit the models
model.fits <- plyr::alply(model.list, 1, function(x,detect.pao){
    cat("\n\n***** Starting ", unlist(x), "\n")
    fit <- RPresence::occMod(model=list(as.formula(paste("psi",x$psi)),
                                        as.formula(paste("p"  ,x$p  )),
                                        as.formula(paste("gamma",x$gamma)),
                                        as.formula(paste("epsilon",x$epsilon))),
          data=detect.pao,type=x$model.type)
    fit <- RPresence.add.derived(fit)
    fit

},detect.pao=cat.pao)


# Look at output from a specified model
model.number <- 3


names(model.fits[[model.number]])
names(model.fits[[model.number]]$real)
model.fits[[model.number]]$beta
names(model.fits[[model.number]]$derived)
model.fits[[model.number]]$derived$psi[1:10,]
model.fits[[model.number]]$real$gamma[1:5,]
model.fits[[model.number]]$real$epsilon[1:5,]

# Estimate of initial occupance
model.fits[[model.number]]$real$psi[grepl('unit1_', row.names(model.fits[[model.number]]$real$psi)),]

# Derived parameters - estimated occupancy for each unit in years 2....
names(model.fits[[model.number]]$derived)
model.fits[[model.number]]$derived$psi[ grepl('unit1_', row.names(model.fits[[model.number]]$derived$psi)),]

# Derived parameters - all of the psi stacked together
model.fits[[model.number]]$derived$all_psi[ grepl('unit1_', row.names(model.fits[[model.number]]$derived$all_psi)),]

# Estimate of  local extinction probability for each unit
model.fits[[model.number]]$real$epsilon[ seq(1, by=nrow(input.history), length.out=length(Nvisits.per.season)-1),]

# Estimate of  local colonization probability for each unit
model.fits[[model.number]]$real$gamma[ seq(1, by=nrow(input.history), length.out=length(Nvisits.per.season)-1),]

# Estimate of probability of detection at each time point for each unit
model.fits[[model.number]]$real$p[ grepl('unit1_', row.names(model.fits[[model.number]]$real$p), fixed=TRUE),]


# Get the change in occupancy
# Not yet possible to estimate the se of these values. May have to use bootstrapping.
model.fits[[model.number]]$derived$lambda [grepl('unit1_', row.names(model.fits[[model.number]]$derived$lambda),  fixed=TRUE),]
model.fits[[model.number]]$derived$lambdap[grepl('unit1_', row.names(model.fits[[model.number]]$derived$lambdap), fixed=TRUE),]



# collect models and make AIC table
aic.table <- RPresence::createAicTable(model.fits)
aic.table$table


# model averaging in the usual way

# initial occupancy
RPresence::modAvg(aic.table, param="psi")[1:5,]

# model averaging of derived parameters such as the occupancy at each time step
ma_all_psi <- RPresence.modAvg.derived(aic.table, param="all_psi")
ma_all_psi[grepl('unit1_', row.names(ma_all_psi),  fixed=TRUE),]


cleanup(ask=FALSE)
