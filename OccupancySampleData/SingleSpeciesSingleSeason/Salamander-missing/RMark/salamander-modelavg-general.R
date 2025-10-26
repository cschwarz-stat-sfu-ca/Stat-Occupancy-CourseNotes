# Blue Ridge Salamander with some added missing values
# Notice how we change the NA to . in the capture history

# Single Species Single Season Occupancy

# Fitting multiple models and model averaging with a general structure


# 2018-08-15 Code contributed by Carl James Schwarz (cschwarz.stat.sfu.cs@gmail.com)

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

library(readxl)
library(RMark)
library(ggplot2)

# get the data read in
# Data for detections should be a data frame with rows corresponding to sites
# and columns to visits.
# The usual 1=detected; 0=not detected; NA=not visited is used.

input.data <- readxl::read_excel("../salamander.xls",
                                 sheet="MissingData",
                                 na='-',
                                 col_names=FALSE)  # notice no column names in row 1 of data file. 

head(input.data)


# Extract the history records and create a capture history
input.history <- data.frame(freq=1,
                            ch=apply(input.data[,1:5],1,paste, collapse=""), stringsAsFactors=FALSE)
head(input.history)

# Change any NA to . in the chapter history
input.history$ch <- gsub("NA",".", input.history$ch, fixed=TRUE)
head(input.history)

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.history)



sal.data <- process.data(data=input.history,
                          model="Occupancy")

summary(sal.data)


# add a survey level covariate. by adding a column to the design matrix

sal.ddl <- make.design.data(sal.data)
sal.ddl

sal.ddl$p$Effort <- factor(c("d1","d1","d2","d2","d2"))
sal.ddl



# What are the parameter names for Single Season Single Species models
setup.parameters("Occupancy", check=TRUE)


# Get the list of models
model.list.csv <- textConnection("
p,            Psi
 ~1,         ~1
 ~time,      ~1
 ~Effort,    ~1")

model.list <- read.csv(model.list.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)
model.list$model.number <- 1:nrow(model.list)
model.list

# fit the models
myobs <- ls()
myobs <- myobs[ grepl("m...",myobs,fixed=TRUE)]
cat("Removing ", myobs, "\n")
rm(list=myobs)

model.fits <- plyr::dlply(model.list, "model.number", function(x,input.data, input.ddl){
  cat("\n\n***** Starting ", unlist(x), "\n")
  #browser()
  
  #fit <- myoccMod(model=list(as.formula(paste("psi",x$psi)),
  fit <- RMark::mark(input.data, ddl=input.ddl,
                     model="Occupancy",
                     model.parameters=list(
                       Psi   =list(formula=as.formula(eval(x$Psi))),
                       p     =list(formula=as.formula(eval(x$p)))
                     )
                     #,brief=TRUE,output=FALSE, delete=TRUE
                     #,invisible=TRUE,output=TRUE  # set for debugging
  )
  mnumber <- paste("m...",formatC(x$model.number, width = 3, format = "d", flag = "0"),sep="")
  assign( mnumber, fit, envir=.GlobalEnv)
  #browser()
  fit
  
},input.data=sal.data, input.ddl=sal.ddl)


# examine individula model results
model.number <-2

summary(model.fits[[model.number]])
model.fits[[model.number]]$results$real
model.fits[[model.number]]$results$beta
model.fits[[model.number]]$results$derived

get.real(model.fits[[model.number]], "Psi", se=TRUE)
get.real(model.fits[[model.number]], "p",    se=TRUE)


# collect models and make AICc table

model.set <- RMark::collect.models( type="Occupancy")
model.set

names(model.set)
model.set$model.table


# model averaged values
get.real(model.set[[1]], "Psi", se=TRUE)
get.real(model.set[[2]], "Psi", se=TRUE)
get.real(model.set[[3]], "Psi", se=TRUE)

Psi.ma <- RMark::model.average(model.set, param="Psi")
Psi.ma


get.real(model.set[[1]], "p", se=TRUE)
get.real(model.set[[2]], "p", se=TRUE)
get.real(model.set[[3]], "p", se=TRUE)

p.ma <- RMark::model.average(model.set, param="p")
p.ma


# cleanup
cleanup(ask=FALSE)

