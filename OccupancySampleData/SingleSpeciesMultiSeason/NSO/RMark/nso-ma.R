# Single Species, Multi Season Occupancy analyais

# Northern Spotted Owl (Strix occidentalis caurina) in California.
# s=55 sites visited up to K=5 times per season between 1997 and 2001 (Y=5).
# Detection probabilities relatively constant within years, but likely different among years.

# Using RMark and several models including different parameterizations of the multi-season model

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

library(plyr)
library(readxl)
library(RMark)
library(ggplot2)

# Get the RMark additional functions 
source(file.path("..","..","..","AdditionalFunctions","RMark.additional.functions.R"))

# get the data read in
# Data for detections should be a data frame with rows corresponding to sites
# and columns to visits.
# The usual 1=detected; 0=not detected; NA=not visited is used.

input.data <- read.csv(file.path("..","NSO.csv"), header=FALSE, skip=2, na.strings="-")
input.data$V1 <- NULL # drop the site number

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.data)
ncol(input.data)
range(input.data, na.rm=TRUE)
sum(is.na(input.data))


# Extract the history records and create a capture history
input.history <- data.frame(freq=1,
                            ch=apply(input.data,1,paste, collapse=""), stringsAsFactors=FALSE)
head(input.history)


# Change any NA to . in the chapter history
input.history$ch <- gsub("NA",".", input.history$ch, fixed=TRUE)
head(input.history)

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.history)


# Create the data structure
max.visit.per.year <- 8
n.season <- 5

nso.data <- process.data(data=input.history, model="RDOccupEG",
                         time.intervals=c( rep( c(rep(0,max.visit.per.year-1),1),n.season-1),
                                           rep(0,max.visit.per.year-1)))

summary(nso.data)

# What are the parameter names for Single Season Single Species models
setup.parameters("RDOccupEG", check=TRUE)

#there are other parameterizations available
setup.parameters("RDOccupPE", check=TRUE) # psi, epsilon, p
setup.parameters("RDOccupPG", check=TRUE) # psi, gamma, p

# time-specific covariates must be added to the ddl's
# look at builtin variables for p
nso.ddl <- make.design.data(nso.data)
nso.ddl$p
str(nso.ddl$p)

nso.ddl$Gamma

# define the list of models to fit
# Notice the commas between the column and the placement of the quotes
# Define the models.
#    model.type can be RDOccupPE, RDOccupPG, RDOccupEG
# The random occupancy model cannot be fit here. See the other code in this directory.
# 
# Notice that the last 3 models are all identical.

model.list.csv <- textConnection("
 p,               Psi,          Gamma,      Epsilon, model.type
~1,              ~1,              ~1,           ~1,       RDOccupEG
~session,        ~1,              ~1,           ~1,       RDOccupEG
~session,        ~1,            ~time,        ~time,      RDOccupEG
~session,        ~time,         ~time,        ~1,         RDOccupPG
~session,        ~time,         ~1,          ~time,       RDOccupPE")


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
  max.visit.per.year <- 8
  n.season <- 5
  
  # we need to process the data in the loop to allow for different data types
  input.data <- process.data(data=input.history, model=x$model.type,
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

names(model.set)
model.set$model.table

# The usual model averaging takes place
# RMark::model.average(model.set, "Psi")  # oops a problem with different model structures

RMark::model.average(model.set, "p") # we can model the average p's

# Model average the derived parameters
# Because RMarks stores psi in different places, we standarize the dervied parameters.
# We need to do this here because collect.models re-extracts the output from MARK and wipes anything else
model.set[-length(model.set)] <- plyr::llply(model.set[-length(model.set)], function(x){RMark.add.derived(x)})
model.set[[1]]$results$derived

RMark.model.average.derived(model.set, "all_psi")

RMark.model.average.derived(model.set, "lambda Rate of Change")

RMark.model.average.derived(model.set, "log odds lambda")

cleanup(ask=FALSE)







