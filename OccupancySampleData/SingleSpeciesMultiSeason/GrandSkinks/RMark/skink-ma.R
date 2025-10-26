# Single Species, MultiSeason Occupancy analyais - multi-model averaging

# Grand Skinks


#   Data has been collected on 352 tors over a 5 year (the seasons) period,^
#   although not all tors (rock piles) were surveyed each year, with up to 3 surveys^
#   of each tor per year.  

#   The 15 columns are in 5 blocks of 3.

#   There is also a site-specific covariate Pasture indicating whether the surrounding^
#   matrix is either predominately the modified habitat (farm pasture, Pasture =1) or
#   "native" grassland (tussock, Pasture = 0).

# 2018-08-18 Code submitted by Carl James Schwarz (cschwarz.stat.sfu.ca@gmail.com)

# Using the RMark package

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

library(car)
library(readxl)
library(RMark)
library(ggplot2)

# Get the RMark additional functions 
source(file.path("..","..","..","AdditionalFunctions","RMark.additional.functions.R"))

# get the data read in
# Data for detections should be a data frame with rows corresponding to sites
# and columns to visits.
# The usual 1=detected; 0=not detected; NA=not visited is used.

input.data <- readxl::read_excel(file.path("..","GrandSkinks.xls"), 
                                    sheet="DetectionHistory",
                                    skip=3, na="-",
                                    col_names=FALSE)
head(input.data)
input.history <- input.data[, 2:16]
head(input.history)

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.history)
ncol(input.history)
range(input.history, na.rm=TRUE)
#Looks like there are some missing values
sum(is.na(input.history))

site.covar <- data.frame(Site=1:nrow(input.data),
                         Pasture=car::recode(input.data$X__18,
                                            "1='Modified'; 0='Native';"))
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
#Five  Seasons, with 3 visits per season
max.visit.per.year <- 3
n.season <- 5

skink.data <- process.data(data=input.history, 
                          model="RDOccupEG",
                          groups = "Pasture",
                          time.intervals=c( rep( c(rep(0,max.visit.per.year-1),1),n.season-1),
                                            rep(0,max.visit.per.year-1)))
summary(skink.data)

# any time specific covariates need to be added to the ddl's in the loop below
# set up the ddls as needed for time-varying covariates.
# you need to do this in the loop because different paraemterizations have different ddl structures


# What are the parameter names for Multi Season Single Species models
setup.parameters("RDOccupEG", check=TRUE)

#there are other parameterizations available
setup.parameters("RDOccupPE", check=TRUE) # psi, epsilon, p
setup.parameters("RDOccupPG", check=TRUE) # psi, gamma, p

# define the list of models to fit
# Notice the commas between the column and the placement of the quotes
# Define the models.
#    model.type can be RDOccupPE, RDOccupPG, RDOccupEG~time
#The random occupancy model cannot be fit here. See the other code in this directory.

model.list.csv <- textConnection("
 p,               Psi,          Gamma,      Epsilon, model.type
~1,              ~1,              ~1,           ~1,       RDOccupEG
~session,           ~1,              ~time,           ~time,       RDOccupEG
~session,           ~Pasture,            ~time,        ~time,      RDOccupEG
~session,           ~Pasture,         ~time*Pasture,           ~time,      RDOccupEG")


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
  n.season <- 5
  
  # we need to process the data in the loop to allow for different data types
  input.data <- process.data(data=input.history,
                             model=x$model.type,
                             time.intervals=c( rep( c(rep(0,max.visit.per.year-1),1),n.season-1),
                                               rep(0,max.visit.per.year-1)),
                             group = "Pasture")
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
model.number <-4
summary(model.fits[[model.number]])
model.fits[[model.number]]$results$real
model.fits[[model.number]]$results$beta
#NOTE: Because a group was used in the process.data step above, the first half
#of each derived parameter table will be for Modified Pastures, and the second
#will be for Native pastures
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
#NOTE: Because a group was used in the process.data step above, the first half
#of each derived parameter table will be for Modified Pastures, and the second
#will be for Native pastures
RMark.model.average.derived(model.set, "all_psi")
RMark.model.average.derived(model.set, "lambda Rate of Change")
RMark.model.average.derived(model.set, "log odds lambda")



# model averaging of derived parameters such as the occupancy at each time step
psi.ma <- RMark.model.average.derived(model.set, "all_psi")
# Note that due to the categorical covariate, there are 10 estimates, but only 5 years
# worth of data. Need to account for this when making Year column
psi.ma$Year <- rep(1:(nrow(psi.ma)/2),2)
psi.ma$parameter <- 'psi'
psi.ma$Pasture = c(rep("Modified",5),rep("Native",5))
psi.ma

# likely more interested in colonization and extinction probabilities
epsilon.ma <- RMark::model.average(model.set, "Epsilon", vcv=TRUE)$estimates
epsilon.ma$Year <- rep(1:(nrow(epsilon.ma)/2),2)
epsilon.ma$parameter <- 'epsilon'
epsilon.ma

gamma.ma <- RMark::model.average(model.set, "Gamma", vcv=TRUE)$estimates
gamma.ma$Year <- rep(1:(nrow(gamma.ma)/2),2)
gamma.ma$parameter <- 'gamma'
gamma.ma

all.est <- plyr::rbind.fill(psi.ma, epsilon.ma, gamma.ma)


ggplot(data=all.est, aes(x=Year,y=estimate, color=parameter))+
  ggtitle("Estimated occupancy, extinction, colonization, over time")+
  geom_point(position=position_dodge(w=0.2))+
  geom_line(position=position_dodge(w=0.2))+
  ylim(0,1)+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.1,position=position_dodge(w=0.2))+
  scale_x_continuous(breaks=1:10)+
  facet_wrap(~Pasture)


cleanup(ask=FALSE)
