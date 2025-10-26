
# Redbreast sunfish occupancy of streams segments over 8 seasons = 4 years * summer*spring
# 3 quadrats were sampled per segment via electrofishing total 24 sampling occasions
# 15 covariates: LN(link magnitude) and standardized Q (streamflow) for the sampling intervals;
# minimum discharge for intervals 1-7, maximum discharge intervals 1-7
  
# Fitting several models using the RMark package
# Single Species Multiple Seasons
# 2018-11-26 code submitted by Neil Faught

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

library(readxl)
library(RMark)
library(ggplot2)

# Get the RMark additional functions 
source(file.path("..","..","..","AdditionalFunctions","RMark.additional.functions.R"))

# get the data read in
# Data for detections should be a data frame with rows corresponding to sites
# and columns to visits.
# The usual 1=detected; 0=not detected; NA=not visited is used.

input.data <- readxl::read_excel(file.path("..","sunfish.xls"), 
                                 skip=1, col_names=TRUE)#No NAs in this data set, make sure to use "na =" statement if there are
head(input.data)

input.history <- input.data[, 1:24]
head(input.history)
# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.history)
ncol(input.history)
range(input.history, na.rm=TRUE)
sum(is.na(input.history))

site.covar <- data.frame(Site=1:nrow(input.data), input.data[,25:ncol(input.data)])
head(site.covar)

#Format the capture history to be used by RMark
input.history <- data.frame(freq=1,
                            ch=apply(input.history,1,paste, collapse=""), stringsAsFactors=FALSE)
head(input.history)

# Change any NA to . in the chapter history
input.history$ch <- gsub("NA",".", input.history$ch, fixed=TRUE)
head(input.history)

#Add site covariates to input history
input.history = cbind(input.history,site.covar)
head(input.history)

#Create the RMark data structure
#8  Seasons, with 3 visits per season

max.visit.per.year <- 3
n.season <- 8

sunfish.data <- process.data(data=input.history, 
                           model="RDOccupEG",
                           time.intervals=c( rep( c(rep(0,max.visit.per.year-1),1),n.season-1),
                                             rep(0,max.visit.per.year-1)))
summary(sunfish.data)

# add visit level covariates to the ddl's in the loop below
# In this case none


# define the list of models to fit
# Notice the commas between the column and the placement of the quotes
# Define the models.
#    model.type can be RDOccupPE, RDOccupPG, RDOccupEG
#The random occupancy model cannot be fit here. See the other code in this directory.

model.list.csv <- textConnection("
                                 p,               Psi,          Gamma,      Epsilon, model.type
                                 ~session,        ~LN,          ~maxQ,      ~minQ,       RDOccupEG
                                 ~session,        ~LN,          ~minQ,      ~minQ,       RDOccupEG
                                 ~session,        ~LN,          ~minQ,      ~maxQ,      RDOccupEG
                                 ~session,        ~1,           ~time,      ~time,      RDOccupEG")


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
  max.visit.per.year <- 3
  n.season <- 8
  
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

#No need to do model averaging in this example as one model has essentially 100% of the weight

# cleanup
cleanup(ask=FALSE)

