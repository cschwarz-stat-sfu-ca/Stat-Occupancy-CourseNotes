# Multi Species, Single Season Occupancy analyais

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

# ***** Only the first "year" of data (2013/2014) prior to any manipulation of the number of cats **** is used

# 2020-01-07 Code contributed by Carl James Schwarz (cschwarz.stat.sfu.cs@gmail.com)

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
# First model occupancy for the rabbits in the first year
# See Table 1 of the paper for the results from a dynamic occupancy

rabbit.detect.vars <- c("Rabbits.2013.2014","X", paste0("X.",1:14))
rabbit.input.history <- input.data[, rabbit.detect.vars]

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(rabbit.input.history)
ncol(rabbit.input.history)
range(rabbit.input.history, na.rm=TRUE)
sum(is.na(rabbit.input.history))
sum(rabbit.input.history, na.rm=TRUE)


other.covar.vars <- c("RabbitPresence_2013_14",  # were rabbits were detected
                      "icats1",    # number of individual cats in year 1
                      "dist2014z", # distance to site with a cat in 2014 
                      "icat2014")  # indicator if cat captured within 500m of the site in 2014
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

# With 16 detections visits, look at the naive occupancy vs. the different covariates
prop.table(xtabs(~RabbitPresence_2013_14+HumanTrail, data=input.history, exclude=NULL, na.action=na.pass), margin=2)
prop.table(xtabs(~RabbitPresence_2013_14+Rabbitat,   data=input.history, exclude=NULL, na.action=na.pass), margin=2)
prop.table(xtabs(~RabbitPresence_2013_14+freshwater, data=input.history, exclude=NULL, na.action=na.pass), margin=2)
prop.table(xtabs(~RabbitPresence_2013_14+salt_butt,  data=input.history, exclude=NULL, na.action=na.pass), margin=2)
prop.table(xtabs(~RabbitPresence_2013_14+icats1,     data=input.history, exclude=NULL, na.action=na.pass), margin=2)

xtabs(~freshwater+salt_butt,   data=input.history, exclude=NULL, na.action=na.pass)



#Create the RMark data structure
# 1 years with 16 camera trapping days/year
max.visit <- 16

rabbit.data <- process.data(data=input.history, 
                                model="Occupancy"
                                )
rabbit.data


# Get the parameter names
# What are the parameter names for Single Season Single Species models
setup.parameters("Occupancy", check=TRUE)


# The top model for the dyanmic occuanpcy is psi(global), .... p(humtrail)
# If you look at the supplemental table for model 20 you see that the "global" terms are
# this is psi(humtrail+rabbitat+devel+fresh+salt_butt+patch+LiDar)

# Notice the commas between the column and the placement of the quotes
# Define the models.
#    model.type can be 
#The random occupancy model cannot be fit here. See the other code in this directory.
# The reserved keywords can be found by looking at the ddl structures earlier

global <- paste0("~",paste0(global.covar.vars, collapse="+"),collapse="")

model.list.csv <- textConnection("
p,               Psi,              model.type, comment
~HumanTrail,   ~global+icats1,       Occupancy, (20) in paper
~HumanTrail,   ~1,                   Occupancy,
~1,            ~global+icats1,       Occupancy,
~1,            ~dist_dev_z,          Occupancy,
~1,            ~Rabbitat,            Occupancy,
~1,            ~freshwater,          Occupancy,
~1,            ~icats1,              Occupancy,
~1,            ~HumanTrail,          Occupancy,
~1,            ~freshwater+icats1,   Occupancy,
~1,            ~freshwater+icats1+Rabbitat,   Occupancy")


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
cleanup(ask=FALSE)

model.fits <- plyr::dlply(model.list, "model.number", function(x,input.history){
  cat("\n\n***** Starting ", unlist(x), "\n")

  # we need to process the data in the loop to allow for different data types
  # notice that max.visits.per.year and n.season are defined outside the function (bad programming practise)
  input.data<- process.data(data=input.history, model=x$model.type)
  # set up the ddls as needed for time-varying covariates.
  # you need to do this in the loop because different paraemterizations have different ddl structures
  input.ddl <- make.design.data(input.data)
  
  model.parameters=list(
    Psi   =list(formula=as.formula(eval(x$Psi))),
    p     =list(formula=as.formula(eval(x$p  )))    )

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

# examine individual model results

model.number <-9

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

# find occupancy for combination of freshwater and icats1
all.ddl <- make.design.data(rabbit.data)
all.ddl$Psi

Temp.df <- expand.grid(freshwater=c(0,1), icats1=0:4)
model.fits[[9]]$results$beta
pred.psi <- covariate.predictions(model.fits[[9]], indices=17, data=Temp.df)
head(pred.psi$estimates)

ggplot(data=pred.psi$estimates, aes(x=icats1, y=estimate, color=as.factor(freshwater)))+
   ggtitle("Estimated rabbit occupancy as a function of cats seen at a site")+
   geom_point()+
   geom_line()+
   geom_ribbon(aes(ymin=lcl, ymax=ucl), alpha=0.1)
  



#####################################################################
#####################################################################
#####################################################################
# Next model static occupancy for the cats
# See Tables S2 and S3 of the paper for dynamic occupancy model

global.covar.vars <- c("HumanTrail","Rabbitat","dist_dev_z","freshwater","salt_butt","hectares_z","LiDar_z")
setdiff(global.covar.vars, names(input.data))


cat.detect.vars <- c("Cats.2013.2014", paste0("X.",31:45))
cat.input.history <- input.data[, cat.detect.vars]

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(cat.input.history)
ncol(cat.input.history)
range(cat.input.history, na.rm=TRUE)
sum(is.na(cat.input.history))
sum(cat.input.history, na.rm=TRUE)


other.covar.vars <- c("Total.CatCaps_2013_14")  # number of cats detected
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


# With 16 detections visits, look at the naive occupancy vs. the different covariates
input.history$CatCaps <- as.numeric(input.history$Total.CatCaps_2013_14>0)

prop.table(xtabs(~CatCaps+HumanTrail, data=input.history, exclude=NULL, na.action=na.pass), margin=2)
prop.table(xtabs(~CatCaps+Rabbitat,   data=input.history, exclude=NULL, na.action=na.pass), margin=2)
prop.table(xtabs(~CatCaps+freshwater, data=input.history, exclude=NULL, na.action=na.pass), margin=2)
prop.table(xtabs(~CatCaps+salt_butt,  data=input.history, exclude=NULL, na.action=na.pass), margin=2)

xtabs(~freshwater+salt_butt,   data=input.history, exclude=NULL, na.action=na.pass)


#Create the RMark data structure

cat.data <- process.data(data=input.history, 
                                model="Occupancy")
cat.data



# all time-specific covariates must be added to the ddl in the loop below


# Get the parameter names
# What are the parameter names for Single Season Single Species models
setup.parameters("Occupancy", check=TRUE)


# The top model is psi(global), gamma(ind change in cats) epsilon (dot) p(humtrail)
# If you look at the supplemental table for model 20 you see that
# this is psi(humtrail+catat+devel+fresh+salt_butt+patch+LiDar)

# Notice the commas between the column and the placement of the quotes

global <- paste0("~",paste0(global.covar.vars, collapse="+"),collapse="")

model.list.csv <- textConnection("
p,     Psi,                       model.type, comment
~1,   ~global,                    Occupancy, (15) Table S2
~1,            ~dist_dev_z,          Occupancy,
~1,            ~Rabbitat,            Occupancy,
~1,            ~freshwater,          Occupancy,
~1,            ~HumanTrail,          Occupancy,
~1,            ~freshwater+Rabbitat, Occupancy,
~1,            ~1,                   Occupancy,
~1,   ~hectares_z+dist_dev_z,                   Occupancy")


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
cleanup(ask=FALSE)

model.fits <- plyr::dlply(model.list, "model.number", function(x,input.history){
  cat("\n\n***** Starting ", unlist(x), "\n")

  # we need to process the data in the loop to allow for different data types
  # notice that max.visits.per.year and n.season are defined outside the function (bad programming practise)
  input.data<- process.data(data=input.history, model=x$model.type)
  # set up the ddls as needed for time-varying covariates.
  # you need to do this in the loop because different paraemterizations have different ddl structures
  input.ddl <- make.design.data(input.data)
  
  model.parameters=list(
    Psi   =list(formula=as.formula(eval(x$Psi))),
    p     =list(formula=as.formula(eval(x$p)))
  )

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


#####################################################################
#####################################################################
#####################################################################
# Finally  model MultiSpecies, Single Season using the 2013/2014 data
# prior to any manipulation of cats

# This was NOT done in the paper and so there is nothing to compare 
# the results with.

# let us look at the raw data of naive occupancy
xtabs(~RabbitPresence+(icats1>0), data=input.data, exclude=NULL, na.action=na.pass)
prop.table(xtabs(~RabbitPresence+(icats1>0), data=input.data, exclude=NULL, na.action=na.pass), margin=1)
prop.table(xtabs(~RabbitPresence+(icats1>0), data=input.data, exclude=NULL, na.action=na.pass), margin=2)


# MARK wants a 2 "digit" code xy for each visit where 
#     x = 0/1 if species A is not detected/detected
# and y = 0/1 if species B is not detected/detected

# It is not clear if we should make the rabbits or the cats
# the "dominant species", so we should try if both ways.
# The 2*species code is the dominant species
input.history2 <- 1*rabbit.input.history + 2*cat.input.history
colnames(input.history2) <- NULL
input.history2[1:2,]


input.chistory <- data.frame(lapply(input.history2, 
                           car::recode, " 0='00'; 1='10';2='01';3='11';", as.numeric=FALSE, as.factor=FALSE),
                           stringsAsFactors=FALSE)
input.history <- data.frame(ch=apply(input.chistory, 1, paste, sep="",collapse=""), freq=1)

# Change any NA to .. in the chapter history
select <- grepl("NA", input.history$ch)
input.history[ select,][1:5,]

input.history$ch <- gsub("NA","..", input.history$ch, fixed=TRUE)
input.history[ select,][1:5,]


other.covar.vars <- NULL  # 
setdiff(other.covar.vars, names(input.data))

site.covar <- input.data[, c(global.covar.vars,other.covar.vars)]
row.names(site.covar) <- NULL
head(site.covar)


#Add site covariates to input history
input.history = cbind(input.history,site.covar)
head(input.history)

catrabbit.data <- process.data(data=input.history,
                         model="2SpecConOccup")  # this parameterization is more stable


# are there any time-specific covariates? If so add them to the ddl's
catrabbit.ddl <- RMark::make.design.data(catrabbit.data)  # need for covariate predictions


# What are the parameter names for Single Season Single Species models
setup.parameters("2SpecConOccup", check=TRUE)


#-----
# Use the psiAB parameterization 
#    Parameters of the model are 
#        psiA  - occupancy probability of species A
#        psiBA - occupancy probability of species B if species A is present
#        psiBa - occupancy probability of species B if species A is absent
#
#    If species are independent thatn psiBA = psiBa.
#       Alternatively, nu = odds(B|A)/odds(B|a) = 1.
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


setup.parameters("2SpecConOccup", check=TRUE)


# Get the list of models. NOtice NO equal signs here
model.list.csv <- textConnection("
PsiA,         PsiBA,           PsiBa,             pA,   pB,   rA,   rBA,   rBa
~1,            ~1,              ~1,               ~1,     ~1,  ~1,   ~1,    ~1
~1,            ~1,              ~1,               ~1,     ~1,  ~1,   ~1,    ~SHARE
~1,            ~1,              ~SHARE,           ~1,     ~1,  ~1,   ~1,    ~1
~1,            ~1,              ~SHARE,           ~1,     ~1,  ~1,   ~1,    ~SHARE")


model.list <- read.csv(model.list.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)
model.list$model.number <- 1:nrow(model.list)
model.list

# fit the models
myobs <- ls()
myobs <- myobs[ grepl("m...",myobs,fixed=TRUE)]
cat("Removing ", myobs, "\n")
rm(list=myobs)
cleanup(ask=FALSE)


model.fits <- plyr::dlply(model.list, "model.number", function(x,input.data, input.ddl){
  cat("\n\n***** Starting ", unlist(x), "\n")
  #browser()
  model.parameters=list(
    PsiA   =list(formula=as.formula(eval(x$PsiA ))),
    PsiBA  =list(formula=as.formula(eval(x$PsiBA))),
    PsiBa  =list(formula=as.formula(eval(x$PsiBa))),
    pA     =list(formula=as.formula(eval(x$pA ))),
    pB     =list(formula=as.formula(eval(x$pB ))),
    rA     =list(formula=as.formula(eval(x$rA ))),
    rBA    =list(formula=as.formula(eval(x$rBA))),
    rBa    =list(formula=as.formula(eval(x$rBa)))
  )
  if(x$rBa == '~SHARE'){
     model.parameters$rBA =list(formula=as.formula(eval(x$rBA)), share=TRUE)
     model.parameters$rBa = NULL
  }
  if(x$PsiBa == '~SHARE'){
     model.parameters$PsiBA =list(formula=as.formula(eval(x$PsiBA)), share=TRUE)
     model.parameters$PsiBa = NULL
  }
 
  fit <- RMark::mark(input.data, ddl=input.ddl,
                     model="2SpecConOccup",
                     model.parameters=model.parameters
                     #,brief=TRUE,output=FALSE, delete=TRUE
                     #,invisible=TRUE,output=TRUE  # set for debugging
          )
  mnumber <- paste("m...",formatC(x$model.number, width = 3, format = "d", flag = "0"),sep="")
  assign( mnumber, fit, envir=.GlobalEnv)
  #browser()
  fit
  
},input.data=catrabbit.data, input.ddl=catrabbit.ddl)


# examine individual model results
model.number <-4

summary(model.fits[[model.number]])
model.fits[[model.number]]$results$real
model.fits[[model.number]]$results$beta
model.fits[[model.number]]$results$derived

get.real(model.fits[[model.number]], "PsiA", se=TRUE)
get.real(model.fits[[model.number]], "PsiBA", se=TRUE)
get.real(model.fits[[model.number]], "PsiBa", se=TRUE)

get.real(model.fits[[model.number]], "pA", se=TRUE)
get.real(model.fits[[model.number]], "pB", se=TRUE)
get.real(model.fits[[model.number]], "rA", se=TRUE)
get.real(model.fits[[model.number]], "rBA", se=TRUE)
get.real(model.fits[[model.number]], "rBa", se=TRUE)

# collect models and make AICc table

model.set <- RMark::collect.models( type="2SpecConOccup")
model.set

names(model.set)
model.set$model.table


# model averaged values at average covariate value

mean(input.data$RabbitPresence)

PsiA.ma <- RMark::model.average(model.set, param="PsiA")
PsiA.ma

PsiBA.ma <- RMark::model.average(model.set, param="PsiBA")
PsiBA.ma

PsiBa.ma <- RMark::model.average(model.set, param="PsiBa")
PsiBa.ma



# look at detection probabiities

RMark::model.average(model.set, param="pA")
RMark::model.average(model.set, param="pB")
RMark::model.average(model.set, param="rA")
RMark::model.average(model.set, param="rBA")
RMark::model.average(model.set, param="rBa")

cleanup(ask=FALSE)
