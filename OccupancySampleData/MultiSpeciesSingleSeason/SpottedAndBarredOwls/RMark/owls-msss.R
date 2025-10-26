# Multi species single season with a list of models

# Using RMark

# Co-occurence of spotted (SO) and barred owls (BO)

# 2018-12-13 code contributed by Neil Faught

#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------

library(car)
library(ggplot2)
library(readxl)
library(reshape2)
library(RMark)

# get the data
SO.data <- readxl::read_excel(file.path("..","SpottedAndBarredOwls.xls"),
                              sheet="RawData", na='-',
                              col_names=FALSE,
                              range = "B11:K161")

BO.data <- readxl::read_excel(file.path("..","SpottedAndBarredOwls.xls"),
                              sheet="RawData", na='-',
                              col_names=FALSE,
                              range = "M11:V161")

# MARK wants a 2 "digit" code xy for each visit where 
#     x = 0/1 if species A is not detected/detected
# and y = 0/1 if species B is not detected/detected

input.data <- 2*BO.data + SO.data
input.data

input.chistory <- data.frame(lapply(input.data, as.character), stringsAsFactors=FALSE)

input.chistory <- data.frame(lapply(input.chistory, 
                                    car::recode, " '0'='00'; '1'='10';'2'='01';'3'='11';", as.numeric=FALSE, as.factor=FALSE),
                             stringsAsFactors=FALSE)
input.history <- data.frame(ch=apply(input.chistory, 1, paste, sep="",collapse=""), freq=1)

# Change any NA to . in the chapter history
select <- grepl("NA", input.history$ch)
input.history[ select,]

input.history$ch <- gsub("NA","..", input.history$ch, fixed=TRUE)
input.history[ select,]

# Read in survey level covariates
night = readxl::read_excel(file.path("..","SpottedAndBarredOwls.xls"),
                           sheet="Nite Covariate", na='-',
                           col_names=FALSE,
                           range = "B11:K161")
head(night)
night[is.na(night)] = 0  # all missing values never used, but RMarks expects a non-missing value when making time varying factor

# recode night as N/D to make your life easier in the long run rather than a simple 0/1 variable
# make each column a factor variable. 
night[] <- lapply(night, car::recode, "0='d'; 1='n'")
night[] <- lapply(night, as.factor)

colnames(night) = paste("night", 1:10, sep="")

# tibbles (from read_excel) make the make.time.factor routine upset. 
night <- data.frame(night)

head(night)
str(night)

# does every night variable have at least two levels?
lapply(night, table)


# paste the time-varying covariate (night) with the input data
night.new <- make.time.factor(night, "night", 1:10, intercept="d",delete=TRUE)

input.history = cbind(input.history, night.new)
head(input.history)

owl.data <- process.data(data=input.history,
                         model="2SpecConOccup")  # this parameterization is more stable

summary(owl.data)

owl.ddl <- RMark::make.design.data(owl.data)  # need for covariate predictions
owl.ddl

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
# Notice that we need to know the levels of night that are created in the input.history from 
# the make.time.factor() function.
model.list.csv <- textConnection("
                                 PsiA,     PsiBA,   PsiBa,    pA,   pB,   rA,   rBA,   rBa
                                 ~1,        ~1,      ~1,     ~nightn,     ~1,  ~1,   ~1, ~SHARE
                                 ~1,        ~1,      ~1,     ~1,     ~1,  ~1,   ~1,    ~1
                                 ~1,        ~1,      ~1,     ~1,     ~1,  ~1,   ~1,    ~SHARE
                                 ~1,        ~1,      ~SHARE, ~1,     ~1,  ~1,   ~1,    ~1
                                 ~1,        ~1,      ~SHARE, ~1,     ~1,  ~1,   ~1,    ~SHARE")



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
  
},input.data=owl.data, input.ddl=owl.ddl)

# examine individual model results
model.number <-5

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

# No need to do model averaging, as 100% of weight is on one model

# cleanup
cleanup(ask=FALSE)

