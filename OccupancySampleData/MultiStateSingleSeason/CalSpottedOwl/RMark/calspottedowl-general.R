# Multi state single season

# Nichols et al (2007) Ecology 88:1395-1400 looked at breeding and
# non-breeding California Spotted Owls in s = 54 sites over k = 5
# surveys.

# Hand-fed mice conrmed occupancy; breeding status only
# confrmed when owl took mouse to nest.
#  . = site not visited.
#  0 = owl not detected.
#  1 = owl detected, but no detection of young.
#  2 = owl detected, along with evidence of young.

# Territories that were expected to be occupied were selectively monitored
# and so psi is expected to be close to 1 and not really of interest.

# See MacKenzie et al, Section 5.7.1 for more details

# 2018-10-01 Code contributed by Carl James Schwarz (cschwarz.stat.sfu.cs@gmail.com)

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
occ.data <- readxl::read_excel(file.path("..","CalSpottedOwl.xls"),
                               sheet="RawData",
                               skip=10, na='.')
head(occ.data)

# create the capture history
input.history <- data.frame(ch  =apply(occ.data[,-1],1,paste,sep="", collapse=""),
                            freq=1)
input.history[1:5,]

# need to deal with missing values
input.history$ch <- gsub("NA",".",input.history$ch)
head(input.history)


# For biological reasons, we believe that errors in identifying states
# varies between the first 2 visits and the last 3 visit (early bs late
# breeding season). We create a site-time covariates.
# and append to the input.history data.frame



#   The covariate matrix would be arranged originally as
#         1 2 3 4 5
#  site1 [E E L L L]
#  site2 [E E L L L]
#   :    [: : : : :]
#  siteN [E E L L L]
#
# and this is stacked by columns


# create the Rmark data file

owl.data <- process.data(data=input.history,
                         model="MSOccupancy")  # this parameterization is more stable

summary(owl.data)

# Time-specific covariates are added to the ddl's.
# We need to modify the design data to create this design matrix
# See the dipper example in Section C.6 of the RMark manual
owl.ddl=make.design.data(owl.data)
names(owl.ddl)
owl.ddl$Delta

owl.ddl$Delta$Per <- factor(c("E","E","L","L","L"))
owl.ddl$Delta
str(owl.ddl$Delta)

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
~1,     ~1,       ~time,   ~1,       ~1
~1,     ~1,       ~1,      ~1,       ~1
~1,     ~1,       ~1,      ~SHARE,   ~Per
~1,     ~1,       ~time,   ~SHARE,   ~Per
~1,     ~1,       ~1,      ~1,       ~Per")

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
  
},input.data=owl.data, input.ddl=owl.ddl)





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


names(model.set)
model.set$model.table


# model averaged values at average covariate value

Psi1.ma <- RMark::model.average(model.set, param="Psi1")
Psi1.ma
Psi2.ma <- RMark::model.average(model.set, param="Psi2")
Psi2.ma

# model average derived parameters
source(file.path("..","..","..","AdditionalFunctions","RMark.additional.functions.R"))
RMark.model.average.derived(model.set, "psi1*psi2")  # Notice that lowercase Psi1 and lowercase Psi2


cleanup(ask=FALSE)
