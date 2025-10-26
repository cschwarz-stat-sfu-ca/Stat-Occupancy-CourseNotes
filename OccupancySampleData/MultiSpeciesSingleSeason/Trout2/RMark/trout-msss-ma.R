# Multi species single season
# Co-occurence of Bull and Brook trout
# Data provided by Parks Canada
# 2018-12-27 code submitted by Neil Faught

library(car)
library(ggplot2)
library(readxl)
library(reshape2)
library(RMark)

# get the data
occ.data <- readxl::read_excel("../Trout2.xlsx",
                               sheet="fish_sp_stacked_juv_delta")
head(occ.data)

# check site code
xtabs(~site+rep, data=occ.data, exclude=NULL, na.action=na.pass)

# some site numbers are "reused", but a combination of watershed x site is unique
nmeasure <- plyr::ddply(occ.data, c("watershed","site"), plyr::summarize, count=length(rep))
nmeasure[ nmeasure$count != 4,]

xtabs(~presence, data=occ.data, exclude=NULL, na.action=na.pass)
occ.data$presence <- as.numeric(occ.data$presence)

# create the detection records for each rep
# Species A = BKTR; Species B=BLTR
# it is believed that BKTR have a preference for warmer water (see later)
# so we code them as the "dominant" species (Species A) so that 
# a linear logisitic regression makes sense for it.
dhist1 <- plyr::ddply(occ.data, c("watershed","site","rep"), plyr::summarize,
                      history=sum(presence*1*(species=="BKTR")+
                                    presence*2*(species=="BLTR")),
                      temp = mean(temp),
                      discharge=mean(discharge))

history <- reshape2::dcast(dhist1, watershed+site+temp+discharge~rep,
                           variable.name="rep",
                           value.var="history")
head(history)
history <- plyr::rename(history, c("1"="r1", "2"="r2"))
xtabs(~r1+r2, data=history, exclude=NULL, na.action=na.pass)

input.data <- history[,c("r1","r2")]
input.data[1:5,]

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

site.covar <- history[,c("watershed","temp","discharge")]
names(site.covar)
head(site.covar)

temp.mean <- mean(site.covar$temp)
temp.sd   <- sd  (site.covar$temp)
site.covar$temp.std <- (site.covar$temp - temp.mean)/ temp.sd

discharge.mean <- mean(site.covar$discharge)
discharge.sd   <- sd  (site.covar$discharge)
site.covar$discharge.std <- (site.covar$discharge - discharge.mean)/ discharge.sd

# Combine site covariates with input history
input.history = cbind(input.history, site.covar)
head(input.history)

# Creat RMark object
trout.data <- process.data(data=input.history,
                           model="2SpecConOccup")  # this parameterization is more stable

summary(trout.data)

# if there are time-specific covariates add them to the ddl's
trout.ddl <- RMark::make.design.data(trout.data)  # need for covariate predictions

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
                                 PsiA,     PsiBA,     PsiBa,       pA,    pB,     rA,   rBA,   rBa
                                 ~1,        ~1,        ~1,         ~1,    ~1,     ~1,   ~1,    ~1
                                 ~1,        ~1,        ~1,         ~1,    ~SHARE, ~1,   ~1,    ~SHARE
                                 ~1,        ~1,        ~SHARE,     ~1,    ~1,     ~1,   ~1,    ~1
                                 ~temp.std, ~temp.std, ~temp.std,  ~1,    ~SHARE, ~1,   ~1, ~SHARE")


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
  if(x$pB == '~SHARE'){
    model.parameters$pA =list(formula=as.formula(eval(x$pA)), share=TRUE)
    model.parameters$pB = NULL
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
  
},input.data=trout.data,input.ddl=trout.ddl)

# examine individual model results
model.number <-4

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

# No need for model averaging as 1 model contains 100% of weight

# cleanup
cleanup(ask=FALSE)
