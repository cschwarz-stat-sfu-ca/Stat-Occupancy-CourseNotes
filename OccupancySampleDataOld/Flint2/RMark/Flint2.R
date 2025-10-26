# Multi species single season with a list of models

# Using RMark

# Flint data set

# 2018-12-13 Code contributed by Neil Faught

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
S1.data <- readxl::read_excel(file.path("..","flint.xlsx"),
                              sheet="Sheet1", na='-',
                              col_names=FALSE,
                              range = "A3:C17")

S2.data <- readxl::read_excel(file.path("..","flint.xlsx"),
                              sheet="Sheet1", na='-',
                              col_names=FALSE,
                              range = "D3:F17")

# MARK wants a 2 "digit" code xy for each visit where 
#     x = 0/1 if species A is not detected/detected
# and y = 0/1 if species B is not detected/detected

input.data <- S1.data + 2*S2.data
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

# Any site covariates?
site.covar <- readxl::read_excel(file.path("..","flint.xlsx"),
                              sheet="Sheet1", na='-',
                              col_names=TRUE,
                              range = "G2:G17")
head(site.covar)
names(site.covar)
str(site.covar)
site.covar$pool = as.factor(site.covar$pool)

input.history = cbind(input.history,site.covar)
head(input.history)

flint.data <- process.data(data=input.history,
                         model="2SpecConOccup",# this parameterization is more stable
                         group = "pool")  

summary(flint.data)

flint.ddl <- RMark::make.design.data(flint.data)  # need for covariate predictions
flint.ddl

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
                                 PsiA,     PsiBA,   PsiBa,    pA,   pB,   rA,   rBA,   rBa
                                 ~1,        ~1,      ~1,     ~1,     ~1,  ~1,   ~1,    ~1
                                 ~1,        ~1,      ~1,     ~1,     ~1,  ~1,   ~1,    ~SHARE
                                 ~1,        ~1,      ~SHARE, ~1,     ~1,  ~1,   ~1,    ~1
                                 ~1,        ~1,      ~SHARE, ~1,     ~1,  ~1,   ~1,    ~SHARE
                                 ~pool,     ~1,      ~SHARE, ~1,     ~1,  ~1,   ~1,    ~SHARE")


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
  
  
  if(x$rBa != '~SHARE' & x$PsiBa != '~SHARE'){
    fit <- RMark::mark(input.data, ddl=input.ddl,
                       model="2SpecConOccup",
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
                       #,brief=TRUE,output=FALSE, delete=TRUE
                       #,invisible=TRUE,output=TRUE  # set for debugging
    )
  }
  if(x$rBa == '~SHARE' & x$PsiBa != '~SHARE'){
    fit <- RMark::mark(input.data, ddl=input.ddl,
                       model="2SpecConOccup",
                       model.parameters=list(
                         PsiA   =list(formula=as.formula(eval(x$PsiA ))),
                         PsiBA  =list(formula=as.formula(eval(x$PsiBA))),
                         PsiBa  =list(formula=as.formula(eval(x$PsiBa))),
                         pA     =list(formula=as.formula(eval(x$pA ))),
                         pB     =list(formula=as.formula(eval(x$pB ))),
                         rA     =list(formula=as.formula(eval(x$rA ))),
                         rBA    =list(formula=as.formula(eval(x$rBA)), share=TRUE)#,
                         #rBa    =list(formula=as.formula(eval(x$rBa)))
                       )
                       #,brief=TRUE,output=FALSE, delete=TRUE
                       #,invisible=TRUE,output=TRUE  # set for debugging
    )
  }
  if(x$rBa != '~SHARE' & x$PsiBa == '~SHARE'){
    fit <- RMark::mark(input.data, ddl=input.ddl,
                       model="2SpecConOccup",
                       model.parameters=list(
                         PsiA   =list(formula=as.formula(eval(x$PsiA ))),
                         PsiBA  =list(formula=as.formula(eval(x$PsiBA)), share=TRUE),
                         #PsiBa  =list(formula=as.formula(eval(x$PsiBa))),
                         pA     =list(formula=as.formula(eval(x$pA ))),
                         pB     =list(formula=as.formula(eval(x$pB ))),
                         rA     =list(formula=as.formula(eval(x$rA ))),
                         rBA    =list(formula=as.formula(eval(x$rBA))),
                         rBa    =list(formula=as.formula(eval(x$rBa)))
                       )
                       #,brief=TRUE,output=FALSE, delete=TRUE
                       #,invisible=TRUE,output=TRUE  # set for debugging
    )
  }
  if(x$rBa == '~SHARE' & x$PsiBa == '~SHARE'){
    fit <- RMark::mark(input.data, ddl=input.ddl,
                       model="2SpecConOccup",
                       model.parameters=list(
                         PsiA   =list(formula=as.formula(eval(x$PsiA ))),
                         PsiBA  =list(formula=as.formula(eval(x$PsiBA)), share=TRUE),
                         #PsiBa  =list(formula=as.formula(eval(x$PsiBa))),
                         pA     =list(formula=as.formula(eval(x$pA ))),
                         pB     =list(formula=as.formula(eval(x$pB ))),
                         rA     =list(formula=as.formula(eval(x$rA ))),
                         rBA    =list(formula=as.formula(eval(x$rBA)), share=TRUE)#,
                         #rBa    =list(formula=as.formula(eval(x$rBa)))
                       )
                       #,brief=TRUE,output=FALSE, delete=TRUE
                       #,invisible=TRUE,output=TRUE  # set for debugging
    )
  }
  mnumber <- paste("m...",formatC(x$model.number, width = 3, format = "d", flag = "0"),sep="")
  assign( mnumber, fit, envir=.GlobalEnv)
  #browser()
  fit
  
},input.data=flint.data, input.ddl = flint.ddl)


# examine individula model results
model.number <-2

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

PsiA.ma <- RMark::model.average(model.set, param="PsiA")
PsiA.ma

# make predictions as a function of pool
# We need to know the indices to the parameters

get.real(model.fits[[model.number]], "PsiA", se=TRUE)
PsiA.ma.e <- RMark::covariate.predictions(model.set, indices=1:4, data=input.history)$estimates
PsiA.ma.e$parameter <- 'PsiA'

get.real(model.fits[[model.number]], "PsiBa", se=TRUE)
PsiBa.ma.e <- RMark::covariate.predictions(model.set, indices=9:12, data=input.history)$estimates
PsiBa.ma.e$parameter <- 'PsiBa'

# Can't go any farther than here due to multiple errors in estimation of regression coefficients.
# Likely there were not enough sites sampled to fit models of this complexity.