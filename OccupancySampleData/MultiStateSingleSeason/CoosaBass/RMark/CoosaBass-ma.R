# Multi state single season
# Coosa Bass collected via electrofishing from four 50-80m sections in streams
# at 54 sites in the Upper Coosa River basin.

# Objective was to evaluate the effect of stream flow variability on Coosa
# Bass reproduction. States are:
#    Species not detected (state = 0)
#    Adult Detected (state = 1)
#    YOY present (state = 2)
# 2 covariates: stream link magnitude (standardized), and CV of stream flow
# during summer.

# Analysis in RMark

# 2018-11-27 code submitted by Neil Faught

# General model averaging framework
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------

library(ggplot2)
library(plyr)
library(readxl)
library(reshape2)
library(RMark)

# get the data
input.data <- readxl::read_excel(file.path("..","CoosaBass.xls"),
                               sheet="Sheet1",
                               skip=13, na='.')
head(input.data)

# create the capture history
input.history <- data.frame(ch  =apply(input.data[,-c(1,6:8)],1,paste,sep="", collapse=""),
                            freq=1)
input.history[1:5,]

# need to deal with missing values
input.history$ch <- gsub("NA",".",input.history$ch)
head(input.history)

# Add site level covariates to input history
site.covar = input.data[,7:8]
nrow(site.covar)
head(site.covar)
names(site.covar) = c("L","CV")
input.history = cbind(input.history,site.covar)
head(input.history)

# create the Rmark data file

bass.data <- process.data(data=input.history,
                         model="MSOccupancy")  # this parameterization is more stable

summary(bass.data)

# time-specific covariates are added to the ddl's
# We need to modify the design data to create this design matrix
# See the dipper example in Section C.6 of the RMark manual
bass.ddl=make.design.data(bass.data)
bass.ddl

setup.parameters("MSOccupancy", check=TRUE)

#-----
# Use the multi-state model, but only for a single season 
# The real parameters are
#    Psi1(i) = probability that a site i is occupied regardless 
#              of reproductive state, Pr(true state = 1 or 2);
#    Psi2(i) = probability that young occurred, given that the site i is occupied,
#              Pr(true state = 2 | true state = 1 or 2);
#
#
#    p1(i, t) = probability that occupancy is detected for site i, 
#               period t, given that true state = 1, Pr(detection | true state = 1);
#    p2(i, t) = probability that occupancy is detected for site i, 
#               period t, given that true state = 2, Pr(detection | true state = 2);
#    Delta(i, t) = probability that evidence of successful reproduction is found, 
#                given detection of occupancy at site i, period t, with successful 
#                reproduction, Pr(classified state 2|true state = 2).

#   Thus, psi1 and psi2 are single parameters relating to a site, 
#   whereas p1, p2, and delta are time-specific parameters, 
#   and thus the number of values in the PIM for each of these parameters equals the number of occasions.
#
#   The product psi1*psi2 is the unconditional probability that a site is 
#   occupied by successful breeders, and is provided as a derived parameter.


# define the list of models to fit
# Notice the commas between the column and the placement of the quotes
# Notice the difference between Time and time in the model terms
model.list.csv <- textConnection("
Psi1,  Psi2,       p1,      p2,    Delta
~1,     ~1,       ~1,      ~SHARE,   ~1
~L,     ~L,       ~L,      ~SHARE,   ~L
~CV,    ~CV,      ~CV,     ~SHARE,   ~CV
")

model.list <- read.csv(model.list.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)
model.list

model.list$model.number <- 1:nrow(model.list)
model.list

# fit the models
myobs <- ls()
myobs <- myobs[ grepl("m...",myobs,fixed=TRUE)]
cat("Removing ", myobs, "\n")
rm(list=myobs)


# fit the models
model.fits <- plyr::dlply(model.list, "model.number", function(x,input.data, input.ddl){
  cat("\n\n***** Starting ", unlist(x), "\n")
  #browser()
  model.parameters=list(
    Psi1   =list(formula=as.formula(eval(x$Psi1 ))),
    Psi2   =list(formula=as.formula(eval(x$Psi2))),
    p1     =list(formula=as.formula(eval(x$p1 ))),
    p2     =list(formula=as.formula(eval(x$p2 ))),
    Delta  =list(formula=as.formula(eval(x$Delta)))
  )
  if(x$p2 == '~SHARE'){
    model.parameters$p1  =list(formula=as.formula(eval(x$p1 )),share=TRUE)
    model.parameters$p2 = NULL
  }
  
  fit <- RMark::mark(input.data, ddl=input.ddl,
                     model="MSOccupancy",
                     model.parameters=model.parameters
                     #,brief=TRUE,output=FALSE, delete=TRUE
                     #,invisible=TRUE,output=TRUE  # set for debugging
  )

  mnumber <- paste("m...",formatC(x$model.number, width = 3, format = "d", flag = "0"),sep="")
  assign( mnumber, fit, envir=.GlobalEnv)
  #browser()
  fit
  
},input.data=bass.data, input.ddl=bass.ddl)



# examine individual model results
model.number <-2

summary(model.fits[[model.number]])
model.fits[[model.number]]$results$real
model.fits[[model.number]]$results$beta
model.fits[[model.number]]$results$derived

get.real(model.fits[[model.number]], "Psi1", se=TRUE)
get.real(model.fits[[model.number]], "Psi2", se=TRUE)

get.real(model.fits[[model.number]], "p1",    se=TRUE)
get.real(model.fits[[model.number]], "p2",    se=TRUE)
get.real(model.fits[[model.number]], "Delta", se=TRUE)

# collect models and make AICc table

model.set <- RMark::collect.models( type="MSOccupancy")
model.set

# No need to do model averaging as nearly 100% of AIC weight is in one model

cleanup(ask=FALSE)
