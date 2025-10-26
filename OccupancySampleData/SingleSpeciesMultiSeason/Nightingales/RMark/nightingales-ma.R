# Single Species, Multi Season Occupancy analyais

#   nightingales over four seasons with covariates on elevation and stream type.
#   We will use only the last 2 years (unsure why)

# 2018-08-15 Code contributed by Carl James Schwarz (cschwarz.stat.sfu.cs@gmail.com)

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

#  RMark package

library(readxl)
library(RMark)
library(ggplot2)

# Get the RMark additional functions 
source(file.path("..","..","..","AdditionalFunctions","RMark.additional.functions.R"))

# get the data read in
# Data for detections should be a data frame with rows corresponding to sites
# and columns to visits.
# The usual 1=detected; 0=not detected; NA=not visited is used.

input.data <- read.csv(file.path("..","nightingales.csv"), header=TRUE, as.is=TRUE, strip.white=TRUE)
head(input.data)

input.history <- input.data[, -1]
head(input.history)


# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.history)
ncol(input.history)
range(input.history, na.rm=TRUE)
sum(is.na(input.history))


# Site covariate. Here there are none
site.covar <- NULL
row.names(site.covar) <- NULL
head(site.covar)

#Format the capture history to be used by RMark
input.history <- data.frame(freq=1,
                            ch=apply(input.history,1,paste, collapse=""), stringsAsFactors=FALSE)
head(input.history)

#Add site covariates to input history
#input.history = cbind(input.history,site.covar)
#head(input.history)

#Create the RMark data structure
# 10 years with 8 visits per year
max.visit.per.year <- 8
n.season <- 10

nightingale.data <- process.data(data=input.history, 
                                model="RDOccupEG",
                                time.intervals=c( rep( c(rep(0,max.visit.per.year-1),1),n.season-1),
                                                  rep(0,max.visit.per.year-1)))
summary(nightingale.data)

# all time-specific covariates must be added to the ddl in the loop below


# Get the parameter names
# What are the parameter names for Single Season Single Species models
setup.parameters("RDOccupEG", check=TRUE)



# define the list of models to fit
# Notice the commas between the column and the placement of the quotes
# Define the models.
#    model.type can be RDOccupPE, RDOccupPG, RDOccupEG~time
#The random occupancy model cannot be fit here. See the other code in this directory.
# The reserved keywords can be found by looking at the ddl structures earlier

model.list.csv <- textConnection("
p,               Psi,          Gamma,      Epsilon,  model.type
~1,              ~1,              ~1,           ~1,       RDOccupEG
~session,        ~1,              ~1,           ~1,       RDOccupEG
~session,        ~1,              ~time,        ~time,    RDOccupEG")


model.list <- read.csv(model.list.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)
model.list$model.number <- 1:nrow(model.list)
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


# examine individula model results

model.number <-2

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

cleanup(ask=FALSE)
