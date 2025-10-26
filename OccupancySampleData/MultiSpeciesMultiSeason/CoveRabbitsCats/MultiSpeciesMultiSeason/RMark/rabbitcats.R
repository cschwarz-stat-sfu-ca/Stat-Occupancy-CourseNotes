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
#
#  Note that the orginal data on the Dryad repository had incorrect cat numbers for 2013/2014 and 2015 data.
#  I wrote to the author and got the corrected data.

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
library(RMark)
library(ggplot2)

# Get the Mark additional functions 
source(file.path("..","..","..","..","AdditionalFunctions","RMark.additional.functions.R"))

# get the data read in
# Data for detections should be a data frame with rows corresponding to sites
# and columns to visits.
# The usual 1=detected; 0=not detected; NA=not visited is used.

input.data <- read.csv(file.path("..","..","Cove_etal_Ecology_Evolution_2018_data-original.csv"), 
                       header=TRUE, as.is=TRUE, strip.white=TRUE, na.string=".")
head(input.data)

# OOPs...# compute the difference in individual cats. These are problematic as the values are the same
# in the two sessions
input.data$CatIndividuals_2013_14
input.data$CatIndividuals_2015
input.data$delta.icats <- input.data$CatIndividuals_2013_14 - input.data$CatIndividuals_2015
input.data$delta.icats

input.data$Total.CatCaps_2013_14
input.data$Total.CatCaps_2015

# here is the updated cat data
input.data2 <- read.csv(file.path("..","..","Cove_etal_Ecology_Evolution_2018_data-updated.csv"), 
                       header=TRUE, as.is=TRUE, strip.white=TRUE, na.string=".")
head(input.data2)
input.data2 <- input.data2[,c("Individualsyear1","Individualsyear2","twoyear_indchange","RabbitPresence")]
input.data <- cbind(input.data, input.data2)  # the second files doesn't have site ids but appears to be sorted in correct order

input.data$CatIndividuals_2013_14 <- NULL
input.data$Total.CatCaps_2015     <- NULL
  
input.data$Individualsyear1
input.data$Individualsyear2
input.data$twoyear_indchange

# Groan, covariate names must be 10 characters or less (a limitation in RMark)
# we need to rename some of the variables
input.data <- plyr::rename(input.data, c("salt_button"           ="salt_butt",
                                         "Individualsyear1"      ="icats1",
                                         "dist_cat2014_z"        ="dist2014z",
                                         "bin_2014cat"           ="icat2014",
                                         "twoyear_indchange"     ="delta_cats"))
names(input.data)

# the paper collapses the different habitats into 3 categories, and creates indicator variables
unique(input.data$DESCRIPTION)
xtabs(~DESCRIPTION+salt_butt,  data=input.data,exclude=NULL, na.action=na.pass)
xtabs(~DESCRIPTION+freshwater, data=input.data,exclude=NULL, na.action=na.pass)
xtabs(~DESCRIPTION+upland,     data=input.data,exclude=NULL, na.action=na.pass)




# distance to sites with cat removed in 2014 vs binary indicator
input.data$icat2014
input.data$dist2014z
ggplot(data=input.data, aes(x=dist2014z, y=icat2014))+
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
                      "dist2014z", # distance to site with a cat in 2014 
                      "icat2014",  # indicator if cat captured within 500m of the site in 2014
                      "delta_cats")# change in number of individual cats captured at a sites
setdiff(other.covar.vars, names(input.data))

site.covar <- input.data[, c(global.covar.vars,other.covar.vars)]
row.names(site.covar) <- NULL
head(site.covar)


#Format the capture history to be used by RMark and deal with missing values
input.history <- data.frame(freq=1,
                            ch=apply(rabbit.input.history,1,paste, collapse=""), stringsAsFactors=FALSE)
head(input.history,n=20)
input.history$ch <- gsub("NA",".",input.history$ch, fixed=TRUE)
head(input.history,n=20)

#Add site covariates to input history
input.history = cbind(input.history,site.covar)
head(input.history)

#Create the RMark data structure
# 2 years with 15 camera trapping days/year
max.visit.per.year <- 16
n.season <- 2

rabbit.data <- process.data(data=input.history, 
                                model="RDOccupEG",
                                time.intervals=c( rep( c(rep(0,max.visit.per.year-1),1),n.season-1),
                                                  rep(0,max.visit.per.year-1)))
rabbit.data


# The paper claims there are 10 observed rabbit colonization events and 7 observed extinction events 
# Is this true? 
dim(rabbit.input.history)
detect.rabbit <- data.frame(detect.r1=apply(rabbit.input.history[, 1:16],1,max, na.rm=TRUE),
                            detect.r2=apply(rabbit.input.history[,17:32],1,max, na.rm=TRUE))
detect.rabbit$ch <- apply(detect.rabbit,1,paste0,collapse="")
head(detect.rabbit)
xtabs(~ch, data=detect.rabbit)


# all time-specific covariates must be added to the ddl in the loop below


# Get the parameter names
# What are the parameter names for Single Season Single Species models
setup.parameters("RDOccupEG", check=TRUE)


# The top model is psi(global), gamma(ind change in cats) epsilon (dot) p(humtrail)
# If you look at the supplemental table for model 20 you see that
# this is psi(humtrail+rabbitat+devel+fresh+salt_butt+patch+LiDar)

# Notice the commas between the column and the placement of the quotes
# Define the models.
#    model.type can be RDOccupPE, RDOccupPG, RDOccupEG~time
#The random occupancy model cannot be fit here. See the other code in this directory.
# The reserved keywords can be found by looking at the ddl structures earlier

global <- paste0("~",paste0(global.covar.vars, collapse="+"),collapse="")

model.list.csv <- textConnection("
p,               Psi,               Gamma,      Epsilon,  model.type, comment
~HumanTrail,   ~global+icats1, ~delta_cats,          ~1,       RDOccupEG, (20) in paper
~HumanTrail,   ~global+icats1,         ~1,           ~1,       RDOccupEG, (13) in paper
~HumanTrail,   ~global+icats1,         ~dist2014z,   ~1,       RDOccupEG
~HumanTrail,   ~global+icats1,         ~icat2014,    ~1,       RDOccupEG
~1,            ~1,         ~1,    ~1,       RDOccupEG")


model.list <- read.csv(model.list.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)
model.list$model.number <- 1:nrow(model.list)
model.list

# substitute for some of the shortcuts
model.list$Psi <- gsub("~global",global, model.list$Psi, )
model.list

# fit the models
myobs <- ls()
myobs <- myobs[ grepl("m...",myobs,fixed=TRUE)]
cat("Removing ", myobs, "\n")
rm(list=myobs)

model.fits <- plyr::dlply(model.list, "model.number", function(x,input.history){
  cat("\n\n***** Starting ", unlist(x), "\n")

  # we need to process the data in the loop to allow for different data types
  # notice that max.visits.per.year and n.season are defined outside the function (bad programming practise)
  input.data<- process.data(data=input.history, model=x$model.type,
                             time.intervals=c( rep( c(rep(0,max.visit.per.year-1),1),n.season-1),
                                               rep(0,max.visit.per.year-1)))
  # set up the ddls as needed for time-varying covariates.
  # you need to do this in the loop because different paraemterizations have different ddl structures
  input.ddl <- make.design.data(input.data)
  
  model.parameters=list(
    Psi   =list(formula=as.formula(eval(x$Psi))),
    p     =list(formula=as.formula(eval(x$p))),
    Epsilon=list(formula=as.formula(eval(x$Epsilon))),
    Gamma  =list(formula=as.formula(eval(x$Gamma)))
  )
  if(x$model.type == "RDOccupPG"){  # psi, gamma, p formulation
    model.parameters$Epsilon= NULL
  }
  if(x$model.type == "RDOccupPE"){  # psi, epsilon, p formulation
    model.parameters$Gamma = NULL
  }  
  
  fit <- RMark::mark(input.data, ddl=input.ddl,
                     model=x$model.type,
                     model.parameters=model.parameters
                     #,brief=TRUE,output=FALSE, delete=TRUE
                     #,invisible=TRUE,output=TRUE  # set for debugging
  )

  mnumber <- paste("m...",formatC(x$model.number, width = 3, format = "d", flag = "0"),sep="")
  assign( mnumber, fit, envir=.GlobalEnv)
  #browser()
  fit
  
},input.history=input.history)
export.MARK(rabbit.data, "Rabbits", model=m...003, replace=TRUE)

# examine individual model results

model.number <-1

summary(model.fits[[model.number]])
model.fits[[model.number]]$results$real
model.fits[[model.number]]$results$beta
model.fits[[model.number]]$results$derived

model.fits[[model.number]]$results$derived$"psi Probability Occupied"
model.fits[[model.number]]$results$derived$"lambda Rate of Change"
model.fits[[model.number]]$results$derived$"log odds lambda"

get.real(model.fits[[model.number]], "Psi", se=TRUE)
get.real(model.fits[[model.number]], "p",    se=TRUE)



# collect models and make AICc table
model.set <- RMark::collect.models( type=NULL)
model.set

# model averaging in the usual way
RMark::model.average(model.set, "Psi")
RMark::model.average(model.set, "p")

# Model average the derived parameters
# Because RMark stores psi in different places, we standarize the dervied parameters.
# We need to do this here because collect.models re-extracts the output from MARK and wipes anything else
model.set[-length(model.set)] <- plyr::llply(model.set[-length(model.set)], function(x){RMark.add.derived(x)})
# NOTE: Because a group was used in the process.data step above, the first half
# of each derived parameter table will be for Prox=No, and the second
# will be for Prox = Yes
model.set[[1]]$results$derived

RMark.model.average.derived(model.set, "all_psi")
RMark.model.average.derived(model.set, "lambda Rate of Change")
RMark.model.average.derived(model.set, "log odds lambda")



# model averaging of derived parameters such as the occupancy at each time step
psi.ma <- RMark.model.average.derived(model.set, "all_psi")
psi.ma$Year <- 1:nrow(psi.ma)
psi.ma$parameter <- 'psi'
psi.ma


# likely more interested in colonization and extinction probabilities
epsilon.ma <- RMark::model.average(model.set, "Epsilon", vcv=TRUE)$estimates
epsilon.ma$Year <- 1:nrow(epsilon.ma)
epsilon.ma$parameter <- 'epsilon'
epsilon.ma

gamma.ma <- RMark::model.average(model.set, "Gamma", vcv=TRUE)$estimates
gamma.ma$Year <- 1:nrow(gamma.ma)
gamma.ma$parameter <- 'gamma'
gamma.ma

all.est <- plyr::rbind.fill(psi.ma, epsilon.ma, gamma.ma)

ggplot(data=all.est, aes(x=Year,y=estimate, color=parameter))+
  ggtitle("Estimated occupancy, extinction, colonization, over time")+
  geom_point(position=position_dodge(w=0.2))+
  geom_line(position=position_dodge(w=0.2))+
  ylim(0,1)+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.1,position=position_dodge(w=0.2))+
  scale_x_continuous(breaks=1:10)




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


#Format the capture history to be used by RMark and deal with missing values
input.history <- data.frame(freq=1,
                            ch=apply(cat.input.history,1,paste, collapse=""), stringsAsFactors=FALSE)
head(input.history,n=20)
input.history$ch <- gsub("NA",".",input.history$ch, fixed=TRUE)
head(input.history,n=20)

#Add site covariates to input history
input.history = cbind(input.history,site.covar)
head(input.history)

#Create the RMark data structure
# 2 years with 15 camera trapping days/year
max.visit.per.year <- 16
n.season <- 2

cat.data <- process.data(data=input.history, 
                                model="RDOccupEG",
                                time.intervals=c( rep( c(rep(0,max.visit.per.year-1),1),n.season-1),
                                                  rep(0,max.visit.per.year-1)))
cat.data


# The paper claims there are 9 observed cat colonization events and 27 observed extinction events 
# Is this true? 
dim(cat.input.history)
detect.cat <- data.frame(detect.r1=apply(cat.input.history[, 1:16],1,max, na.rm=TRUE),
                            detect.r2=apply(cat.input.history[,17:32],1,max, na.rm=TRUE))
detect.cat$ch <- apply(detect.cat,1,paste0,collapse="")
head(detect.cat)
xtabs(~ch, data=detect.cat)


# all time-specific covariates must be added to the ddl in the loop below


# Get the parameter names
# What are the parameter names for Single Season Single Species models
setup.parameters("RDOccupEG", check=TRUE)


# The top model is psi(global), gamma(ind change in cats) epsilon (dot) p(humtrail)
# If you look at the supplemental table for model 20 you see that
# this is psi(humtrail+catat+devel+fresh+salt_butt+patch+LiDar)

# Notice the commas between the column and the placement of the quotes
# Define the models.
#    model.type can be RDOccupPE, RDOccupPG, RDOccupEG~time
#The random occupancy model cannot be fit here. See the other code in this directory.
# The reserved keywords can be found by looking at the ddl structures earlier

global <- paste0("~",paste0(global.covar.vars, collapse="+"),collapse="")

model.list.csv <- textConnection("
p,                                 Psi,               Gamma,      Epsilon,  model.type, comment
~HumanTrail+Rabbitat+year_2013,   ~global,                ~1,      ~year_2013,       RDOccupEG, (15) Table S2")


model.list <- read.csv(model.list.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)
model.list$model.number <- 1:nrow(model.list)
model.list

# substitute for some of the shortcuts
model.list$Psi <- gsub("~global",global, model.list$Psi, )
model.list

# fit the models
myobs <- ls()
myobs <- myobs[ grepl("m...",myobs,fixed=TRUE)]
cat("Removing ", myobs, "\n")
rm(list=myobs)

model.fits <- plyr::dlply(model.list, "model.number", function(x,input.history){
  cat("\n\n***** Starting ", unlist(x), "\n")

  # we need to process the data in the loop to allow for different data types
  # notice that max.visits.per.year and n.season are defined outside the function (bad programming practise)
  input.data<- process.data(data=input.history, model=x$model.type,
                             time.intervals=c( rep( c(rep(0,max.visit.per.year-1),1),n.season-1),
                                               rep(0,max.visit.per.year-1)))
  # set up the ddls as needed for time-varying covariates.
  # you need to do this in the loop because different paraemterizations have different ddl structures
  input.ddl <- make.design.data(input.data)
  
  model.parameters=list(
    Psi   =list(formula=as.formula(eval(x$Psi))),
    p     =list(formula=as.formula(eval(x$p))),
    Epsilon=list(formula=as.formula(eval(x$Epsilon))),
    Gamma  =list(formula=as.formula(eval(x$Gamma)))
  )
  if(x$model.type == "RDOccupPG"){  # psi, gamma, p formulation
    model.parameters$Epsilon= NULL
  }
  if(x$model.type == "RDOccupPE"){  # psi, epsilon, p formulation
    model.parameters$Gamma = NULL
  }  
  
  fit <- RMark::mark(input.data, ddl=input.ddl,
                     model=x$model.type,
                     model.parameters=model.parameters,
                     se=TRUE
                     #,brief=TRUE,output=FALSE, delete=TRUE
                     #,invisible=TRUE,output=TRUE  # set for debugging
  )

  mnumber <- paste("m...",formatC(x$model.number, width = 3, format = "d", flag = "0"),sep="")
  assign( mnumber, fit, envir=.GlobalEnv)
  #browser()
  fit
  
},input.history=input.history)
export.MARK(cat.data, "cats", model=m...001, replace=TRUE)

# examine individula model results

model.number <-1

summary(model.fits[[model.number]])
model.fits[[model.number]]$results$real
model.fits[[model.number]]$results$beta
model.fits[[model.number]]$results$derived

model.fits[[model.number]]$results$derived$"psi Probability Occupied"
model.fits[[model.number]]$results$derived$"lambda Rate of Change"
model.fits[[model.number]]$results$derived$"log odds lambda"

get.real(model.fits[[model.number]], "Psi", se=TRUE)
get.real(model.fits[[model.number]], "p",    se=TRUE)



# collect models and make AICc table
model.set <- RMark::collect.models( type=NULL)
model.set

# model averaging in the usual way
RMark::model.average(model.set, "Psi")
RMark::model.average(model.set, "p")

# Model average the derived parameters
# Because RMark stores psi in different places, we standarize the dervied parameters.
# We need to do this here because collect.models re-extracts the output from MARK and wipes anything else
model.set[-length(model.set)] <- plyr::llply(model.set[-length(model.set)], function(x){RMark.add.derived(x)})
model.set[[1]]$results$derived

RMark.model.average.derived(model.set, "all_psi")
RMark.model.average.derived(model.set, "lambda Rate of Change")
RMark.model.average.derived(model.set, "log odds lambda")


#####################################################################
#####################################################################
#####################################################################
# Finally  model dynamic occupancy for both species
# See Table 3 in the paper

# MARK wants a 2 "digit" code xy for each visit where 
#     x = 0/1 if species A is not detected/detected
# and y = 0/1 if species B is not detected/detected

input.history2 <- 2*rabbit.input.history + 1*cat.input.history
colnames(input.history2) <- NULL
input.history2[1:2,]

input.chistory <- data.frame(lapply(as.data.frame(input.history2), as.character), stringsAsFactors=FALSE)

input.chistory <- data.frame(lapply(input.chistory, 
                           car::recode, " '0'='00'; '1'='10';'2'='01';'3'='11';", as.numeric=FALSE, as.factor=FALSE),
                           stringsAsFactors=FALSE)
input.history <- data.frame(ch=apply(input.chistory, 1, paste, sep="",collapse=""), freq=1)

# Change any NA to .. in the chapter history
select <- grepl("NA", input.history$ch)
input.history[ select,][1:5,]

input.history$ch <- gsub("NA","..", input.history$ch, fixed=TRUE)
input.history[ select,][1:5,]


other.covar.vars <- c("year_2013")  # indicator if site first measured in 2013 or 2014
setdiff(other.covar.vars, names(input.data))

site.covar <- input.data[, c(global.covar.vars,other.covar.vars)]
row.names(site.covar) <- NULL
head(site.covar)


#Add site covariates to input history
input.history = cbind(input.history,site.covar)
head(input.history)

#Create the RMark data structure
# 2 years with 15 camera trapping days/year
max.visit.per.year <- 16
n.season <- 2

both.data <- process.data(data=input.history, 
                                model="RD2SpGEConOcc",
                                time.intervals=c( rep( c(rep(0,max.visit.per.year-1),1),n.season-1),
                                                  rep(0,max.visit.per.year-1)))
both.data


# What are the parameter names for Single Season Single Species models
setup.parameters("RD2SpGEConOcc", check=TRUE)


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
# The problem is that RMark doesn't allow you share the GammaA and GammaAB terms etc.
# So we need to construct a base model to export to MARK and do the bulk of the analysis there.
# Here A=cat, B=rabbit
model.list.csv <- textConnection("
PsiA,  PsiBA,   PsiBa,   GammaA, GammaB, GammaAB, GammaBA,   EpsilonA, EpsilonB, EpsilonAB, EpsilonBA,    pA, pB, rA, rBA,rBa
~dh,    ~fw+ps, ~fw+ps, ~y2013,  ~1,     ~y2013,  ~1,        ~y2013,    ~1,      ~y2013,    ~1,            ~1, ~1, ~1, ~1, ~1
~dh,    ~fw+ps, ~SHARE, ~1,      ~1,     ~1,      ~1,        ~1,       ~1,      ~1,         ~1,            ~1, ~1, ~1, ~1, ~1
~1,      ~1,    ~1,     ~1,      ~1,     ~1,      ~1,         ~1,       ~1,      ~1,         ~1,           ~1, ~1, ~1, ~1, ~SHARE
~UR,     ~UR,   ~UR,   ~UR,      ~UR,    ~UR,   ~UR,        ~UR,       ~UR,    ~UR,         ~UR,          ~1, ~1, ~1, ~1, ~1")


model.list <- read.csv(model.list.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)
model.list$model.number <- 1:nrow(model.list)
model.list <- model.list[1:3,]

model.list$PsiA   <- gsub("fw","freshwater",model.list$PsiA )
model.list$PsiBA  <- gsub("fw","freshwater",model.list$PsiBA)
model.list$PsiBa  <- gsub("fw","freshwater",model.list$PsiBa)

model.list$PsiA   <- gsub("ps","hectares_z",model.list$PsiA )
model.list$PsiBA  <- gsub("ps","hectares_z",model.list$PsiBA)  # patch size
model.list$PsiBa  <- gsub("ps","hectares_z",model.list$PsiBa)  # patch size

model.list$GammaA  <- gsub("y2013","year_2013", model.list$GammaA )
model.list$GammaAB <- gsub("y2013","year_2013", model.list$GammaAB)
model.list$PsiA   <- gsub("dh","dist_dev_z",model.list$PsiA)  # distance to development



# are there any time varying covariates that need to be added to the ddls? 
both.ddl <- RMark::make.design.data(both.data)  # need for covariate predictions


# fit the models
myobs <- ls()
myobs <- myobs[ grepl("m...",myobs,fixed=TRUE)]
cat("Removing ", myobs, "\n")
rm(list=myobs)

model.fits <- plyr::dlply(model.list[1,,drop=FALSE], "model.number", function(x,input.data, input.ddl){
  cat("\n\n***** Starting ", unlist(x), "\n")
  #browser()
  
  # Now we need to set up the model parameters accounting for sharing
  model.parameters= list(
    PsiA     =list(formula=as.formula(eval(x$PsiA ))),
    PsiBA    =list(formula=as.formula(eval(x$PsiBA))),
    PsiBa    =list(formula=as.formula(eval(x$PsiBa))),
    GammaA   =list(formula=as.formula(eval(x$GammaA))),
    GammaB   =list(formula=as.formula(eval(x$GammaB))),
    GammaAB  =list(formula=as.formula(eval(x$GammaAB))),
    GammaBA  =list(formula=as.formula(eval(x$GammaBA))),
    EpsilonA =list(formula=as.formula(eval(x$EpsilonA))),
    EpsilonB =list(formula=as.formula(eval(x$EpsilonB))),
    EpsilonAB=list(formula=as.formula(eval(x$EpsilonAB))),
    EpsilonBA=list(formula=as.formula(eval(x$EpsilonBA))),
    pA       =list(formula=as.formula(eval(x$pA ))),
    pB       =list(formula=as.formula(eval(x$pB ))),
    rA       =list(formula=as.formula(eval(x$rA ))),
    rBA      =list(formula=as.formula(eval(x$rBA))),
    rBa      =list(formula=as.formula(eval(x$rBa)))
  )
  if(x$PsiBa == "~SHARE"){
    model.parameters$PsiBA  =list(formula=as.formula(eval(x$PsiBA)), share=TRUE)
    model.parameters$PsiBa  =NULL
  }
  #if(x$GammaA == "~SHARE"){
  #  model.parameters$GammaAB  =list(formula=as.formula(eval(x$GammaAB)), share=TRUE)
  #  model.parameters$GammaA  =NULL
  #}
  if(x$GammaA=="~SHARE"){
    model.parameters$EpsilonA  =list(formula=as.formula(eval(x$EpsilonA)), share=TRUE)
    model.parameters$GammaA  =NULL
  }
  if(x$rBa == '~SHARE'){
    model.parameters$rBA = list(formula=as.formula(eval(x$rBA)), share=TRUE)
    model.parameters$rBa = NULL
  }
  print(model.parameters)
  fit <- RMark::mark(input.data, ddl=input.ddl,
                     model="RD2SpGEConOcc",
                     model.parameters=model.parameters
                     #,brief=TRUE,output=FALSE, delete=TRUE
                     #,invisible=TRUE,output=TRUE  # set for debugging
  )
  mnumber <- paste("m...",formatC(x$model.number, width = 3, format = "d", flag = "0"),sep="")
  assign( mnumber, fit, envir=.GlobalEnv)
  #browser()
  fit
  
},input.data=both.data, input.ddl=both.ddl)
export.MARK(both.data, "both", model=m...001, replace=TRUE)


# examine individual model results

model.number <-1

summary(model.fits[[model.number]])
model.fits[[model.number]]$results$real
model.fits[[model.number]]$results$beta
model.fits[[model.number]]$results$derived

model.fits[[model.number]]$results$derived$"psi Probability Occupied"
model.fits[[model.number]]$results$derived$"lambda Rate of Change"
model.fits[[model.number]]$results$derived$"log odds lambda"

get.real(model.fits[[model.number]], "Psi", se=TRUE)
get.real(model.fits[[model.number]], "p",    se=TRUE)



# collect models and make AICc table
model.set <- RMark::collect.models( type=NULL)
model.set



# make predictions for rabbit and cat initial occupancy
# We need to know the all.diff.index indices to the parameters
# Because we don't have any time effects in the model, we only need to find predictions for one session
setup.parameters("RD2SpGEConOcc", check=TRUE)

model.number <- 2

# initial occupancy of cats
get.real(model.fits[[model.number]], "PsiA", se=TRUE)
PsiA <- RMark::covariate.predictions(model.set, indices=1, data=input.history)$estimates
PsiA$parameter <- "PsiA"

plotdata <- plyr::rbind.fill(PsiA)
ggplot(data=plotdata, aes(x=dist_dev_z, y=estimate))+
   ggtitle("Effects of distance to development for initial occupancy of cats")+
   geom_point()+
   geom_ribbon(aes(ymin=lcl, ymax=ucl), alpha=0.2)+
   facet_wrap(~parameter, ncol=2)



get.real(model.fits[[model.number]], "PsiBA", se=TRUE)
PsiBA <- RMark::covariate.predictions(model.fits[[model.number]], indices=2, data=input.history)$estimates
PsiBA$parameter <- "PsiBA"

get.real(model.fits[[model.number]], "PsiBa", se=TRUE)
PsiBa <- RMark::covariate.predictions(model.fits[[model.number]], indices=3, data=input.history, drop=FALSE)$estimates
PsiBa$parameter <- "PsiBa"

head(PsiBA)
head(PsiBa)

plotdata <- plyr::rbind.fill(PsiBA, PsiBa)
ggplot(data=plotdata, aes(x=hectares_z, y=estimate, color=as.factor(freshwater)))+
   ggtitle("Effects of patch size for initial occupancy of rabbits (B) \nwith/without cats(A/a)")+
   geom_point()+
   geom_ribbon(aes(ymin=lcl, ymax=ucl), alpha=0.2)+
   facet_wrap(~parameter, ncol=2)




get.real(model.fits[[model.number]], "GammaA", se=TRUE)
GammaA <- RMark::covariate.predictions(model.set, indices=4, data=input.history)$estimates
GammaA$parameter <- 'GammaA'

get.real(model.fits[[model.number]], "GammaB", se=TRUE)
GammaB <- RMark::covariate.predictions(model.set, indices=16, data=input.history)$estimates
GammaB$parameter <- 'GammaB'

get.real(model.fits[[model.number]], "GammaAB", se=TRUE)
GammaAB <- RMark::covariate.predictions(model.set, indices=28, data=input.history)$estimates
GammaAB$parameter <- 'GammaAB'

get.real(model.fits[[model.number]], "GammaBA", se=TRUE)
GammaBA <- RMark::covariate.predictions(model.set, indices=40, data=input.history)$estimates
GammaBA$parameter <- 'GammaBA'


plotdata <- plyr::rbind.fill(GammaA, GammaB, GammaAB, GammaBA)
ggplot(data=plotdata, aes(x=UR, y=estimate))+
   ggtitle("Effects of urbanization on local colonization")+
   geom_point()+
   geom_ribbon(aes(ymin=lcl, ymax=ucl), alpha=0.2)+
   facet_wrap(~parameter, ncol=2)



# look at the extinction probabilities
get.real(model.fits[[model.number]], "EpsilonA", se=TRUE)
EpsilonA <- RMark::covariate.predictions(model.set, indices=52, data=input.history)$estimates
EpsilonA$parameter <- 'EpsilonA'

get.real(model.fits[[model.number]], "EpsilonB", se=TRUE)
EpsilonB <- RMark::covariate.predictions(model.set, indices=64, data=input.history)$estimates
EpsilonB$parameter <- 'EpsilonB'

get.real(model.fits[[model.number]], "EpsilonAB", se=TRUE)
EpsilonAB <- RMark::covariate.predictions(model.set, indices=76, data=input.history)$estimates
EpsilonAB$parameter <- 'EpsilonAB'

get.real(model.fits[[model.number]], "EpsilonBA", se=TRUE)
EpsilonBA <- RMark::covariate.predictions(model.set, indices=88, data=input.history)$estimates
EpsilonBA$parameter <- 'EpsilonBA'


plotdata <- plyr::rbind.fill(EpsilonA, EpsilonB, EpsilonAB, EpsilonBA)

ggplot(data=plotdata, aes(x=UR, y=estimate))+
  ggtitle("Effects of urbanization on local extinction")+
  geom_point()+
  geom_ribbon(aes(ymin=lcl, ymax=ucl), alpha=0.2)+
  facet_wrap(~parameter, ncol=2)






cleanup(ask=FALSE)
