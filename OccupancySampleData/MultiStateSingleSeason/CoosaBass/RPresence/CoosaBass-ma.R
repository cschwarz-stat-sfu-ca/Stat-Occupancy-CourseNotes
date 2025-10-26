# Multi state single season
# Coosa Bass collected via electrofishing from four 50-80m sections in streams
# at 54 sites in the Upper Coosa River basin. 
# Objective was to evaluate the effect of stream flow variability on Coosa
# Bass reproduction. States are:
# Species not detected (state = 0)
# Adult Detected (state = 1)
# YOY present (state = 2)
# 2 covariates: stream link magnitude (standardized), and CV of stream flow
# during summer.

# Analysis in RPresence

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
library(RPresence)

# get the data
input.data <- readxl::read_excel(file.path("..","CoosaBass.xls"),
                                 sheet="Sheet1",
                                 skip=13, na='.')
head(input.data)


input.history <- input.data[,2:5]
input.history[1:5,]
nrow(input.history)

Nsites <- nrow(input.history)
Nvisits<- ncol(input.history)

# Extract site level covariates from input data
unit.cov = input.data[,7:8]
names(unit.cov) = c("L","CV")
head(unit.cov)

# create the pao file
bass.pao <- createPao(input.history,
                    unitcov = unit.cov,
                    title="bass multi state single season")
summary(bass.pao)

#-----
# Use the multi-state model, but only for a single season 

# define the list of models to fit
# Notice the commas between the column and the placement of the quotes


#  psi	occupancy probability,
# r	 probability of being in the second state, conditional on the unit being occupied. 
# p	 detection probability. .
# delta	probability of detecting the second state in a survey, conditional on the species being detected in the survey. 
model.list.csv <- textConnection("
psi,     r,        p,     delta
~1,     ~1,       ~1,      ~1
~L,     ~L,       ~L,      ~L
~CV,    ~CV,      ~CV,     ~CV")

model.list <- read.csv(model.list.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)
model.list


# fit the models
model.fits <- plyr::alply(model.list, 1, function(x,detect.pao){
  cat("\n\n***** Starting ", unlist(x), "\n")
  fit <- RPresence::occMod(model=list(as.formula(paste("psi",  x$psi)),
                                      as.formula(paste("r"    ,x$r  )),
                                      as.formula(paste("p"    ,x$p  )),
                                      as.formula(paste("delta",x$delta))),
                           data=detect.pao,type="do.ms.2")
  fit
},detect.pao=bass.pao)




# Look the output from a specific model
check.model <- 3

names(model.fits[[check.model]])
model.fits[[check.model]]$beta

names(model.fits[[check.model]]$real)
model.fits[[check.model]]$real$psi[1:5,]
model.fits[[check.model]]$real$Cpsi0[1:5,]
model.fits[[check.model]]$real$Cpsi1[1:5,]
model.fits[[check.model]]$real$Cpsi2[1:5,]
model.fits[[check.model]]$real$r[1:5,]
model.fits[[check.model]]$real$CR0[1:5,]
model.fits[[check.model]]$real$CR1[1:5,]
model.fits[[check.model]]$real$CR2[1:5,]
model.fits[[check.model]]$real$delta[1:5,]

model.fits[[check.model]]$real$p1[1:5,]
model.fits[[check.model]]$real$p2[1:5,]

names(model.fits[[check.model]]$derived)


# Model averaging
aic.table <- RPresence::createAicTable(model.fits)
aic.table$table

names(aic.table)

# No need to do model averaging as nearly 100% of AIC weight is in one model