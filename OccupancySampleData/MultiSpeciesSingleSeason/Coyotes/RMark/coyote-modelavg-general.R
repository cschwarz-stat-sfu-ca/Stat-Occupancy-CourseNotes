# Multi species single season with a list of models

# Using RMark

# Camera trap detection / non-detection data for bobcats, coyotes, gray foxes, 
# and red foxes from 6 mid-Atlantic states and the District of Columbia, USA; 
# associated site-level covariates.

# Rota CT, Ferreira MAR, Kays RW, Forrester TD, Kalies EL, McShea WJ, Parsons AW, Millspaugh JJ (2016) 
# A multispecies occupancy model for two or more interacting species.
# Methods in Ecology and Evolution 7(10): 1164-1173. 
# https://doi.org/10.1111/2041-210x.12587
#
# We will only look at selected pairs of species

# 2018-12-15e contributed by Carl James Schwarz (cschwarz.stat.sfu.cs@gmail.com)


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
# importing detection / non-detection data for each of the 4 species at each camera trap.
bobcat <- as.matrix(read.csv(file.path("..","Data",'Bobcat.csv'  ), header=FALSE, as.is=TRUE, strip.white=TRUE))
coyote <- as.matrix(read.csv(file.path("..","Data",'Coyote.csv'  ), header=FALSE, as.is=TRUE, strip.white=TRUE))
gryfox <- as.matrix(read.csv(file.path("..","Data",'Gray Fox.csv'), header=FALSE, as.is=TRUE, strip.white=TRUE))
redfox <- as.matrix(read.csv(file.path("..","Data",'Red Fox.csv' ), header=FALSE, as.is=TRUE, strip.white=TRUE))


# MARK wants a 2 "digit" code xy for each visit where 
#     x = 0/1 if species A is not detected/detected
# and y = 0/1 if species B is not detected/detected

input.data <- 2*coyote + 1*bobcat
input.data

input.chistory <- data.frame(lapply(as.data.frame(input.data), as.character), stringsAsFactors=FALSE)

input.chistory <- data.frame(lapply(input.chistory, 
                           car::recode, " '0'='00'; '1'='10';'2'='01';'3'='11';", as.numeric=FALSE, as.factor=FALSE),
                           stringsAsFactors=FALSE)
input.history <- data.frame(ch=apply(input.chistory, 1, paste, sep="",collapse=""), freq=1)


# Change any NA to . in the chapter history
select <- grepl("NA", input.history$ch)
input.history[ select,][1:5,]

input.history$ch <- gsub("NA","..", input.history$ch, fixed=TRUE)
input.history[ select,][1:5,]


# total number of camera days
cday <- apply(bobcat, 1, function(x) sum(is.na(x) == F))


# Any site covariates?
site.covar <- read.csv(file.path("..","Data","psi covariates.csv"), header=TRUE, as.is=TRUE, strip.white=TRUE)
names(site.covar)
names(site.covar) <- make.names(names(site.covar))
names(site.covar)


# Standardize site.covariates

site.covar$Dist_std        = scale(site.covar$Dist_5km)
site.covar$HDens_5km_std   = scale(site.covar$HDens_5km)
site.covar$Latitute_std    = scale(site.covar$Latitude)
site.covar$Longitude_std   = scale(site.covar$Longitude)
site.covar$People_site_std = scale(site.covar$People_site)

# There is also a covariate for each site about the maximum distance a camera can see
site.covar2 <- read.csv(file.path("..","Data","p covariates.csv"), header=TRUE, as.is=TRUE, strip.white=TRUE)
names(site.covar2)
names(site.covar2) <- make.names(names(site.covar2))
names(site.covar2)

site.covar2$Det_std <- scale(site.covar2$Det_dist)

# add the site covariate to the chapture hisotry
input.history <- cbind(input.history, site.covar, site.covar2)
head(input.history)

Nsites <- nrow(input.history)
Nvisits<- nchar(input.history$ch[1])


sal.data <- process.data(data=input.history,
                         model="2SpecConOccup")  # this parameterization is more stable

summary(sal.data)

# are there any time-specific covariates? If so add them to the ddl's
sal.ddl <- RMark::make.design.data(sal.data)  # need for covariate predictions


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
PsiA,         PsiBA,           PsiBa,             pA,   pB,   rA,   rBA,   rBa
~Dist_std,      ~Dist_std,      ~Dist_std,        ~1,     ~1,  ~1,   ~1,    ~SHARE
~1,            ~1,              ~1,               ~1,     ~1,  ~1,   ~1,    ~1
~1,            ~1,              ~1,               ~1,     ~1,  ~1,   ~1,    ~SHARE
~1,            ~1,              ~SHARE,           ~1,     ~1,  ~1,   ~1,    ~1
~1,            ~1,              ~SHARE,           ~1,     ~1,  ~1,   ~1,    ~SHARE")


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
  
},input.data=sal.data, input.ddl=sal.ddl)


# examine individual model results
model.number <-1

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

# make predictions as a function of distance
# We need to know the indices to the parameters
# Look at the all.diff.index values to get thei correct index number


get.real(model.fits[[model.number]], "PsiA", se=TRUE)
PsiA.ma.p <- RMark::covariate.predictions(model.set, indices=1, data=input.history)$estimates
PsiA.ma.p$parameter <- 'PsiA'

get.real(model.fits[[model.number]], "PsiBa", se=TRUE)
PsiBa.ma.p <- RMark::covariate.predictions(model.set, indices=3, data=input.history)$estimates
PsiBa.ma.p$parameter <- 'PsiBa'

get.real(model.fits[[model.number]], "PsiBA", se=TRUE)
PsiBA.ma.p <- RMark::covariate.predictions(model.set, indices=2, data=input.history)$estimates
PsiBA.ma.p$parameter <- 'PsiBA'

# Not easy to get the se of dervied parameters at the covariate values
PsiB.ma.p <- data.frame(estimate=PsiBA.ma.p$estimate*PsiA.ma.p$estimate +
                                 PsiBa.ma.p$estimate*(1-PsiA.ma.p$estimate),
                        Dist_5km=PsiA.ma.p$Dist_5km, stringsAsFactors=FALSE)
PsiB.ma.p$parameter <- "PsiB"


plotdata <- plyr::rbind.fill(PsiA.ma.p,
                             PsiBA.ma.p,
                             PsiBa.ma.p,
                             PsiB.ma.p)

# It is clear that the fit is not satisfactory with a huge se for PsiBa and very steep decline.
# Try reversing the species and you get a much better fit
ggplot(data=plotdata, aes(x=Dist_5km, y=estimate))+
   ggtitle("Effects of proportion of disturbance within 5 km")+
   geom_point()+
   geom_ribbon(aes(ymin=lcl, ymax=ucl), alpha=0.2)+
   facet_wrap(~parameter, ncol=2)



# Why is psiBa so bad? This may be a case of complete separation.
# The detection probabilities are all fairly large so that over 5 visits, there is almost certain detection
# if the species is present. Consequently, the observed occupancy is likely a pretty good indication of the 
# actual occupancy.
RMark::model.average(model.set, param="pA")
RMark::model.average(model.set, param="pB")
RMark::model.average(model.set, param="rA")
RMark::model.average(model.set, param="rBA")
RMark::model.average(model.set, param="rBa")

# Compute the observed occupance of each species
obs.occ <- data.frame(Dist_5km = site.covar$Dist_5km,
                      ooC = apply(coyote, 1,max, na.rm=TRUE),
                      ooB = apply(bobcat ,1,max, na.rm=TRUE))
head(obs.occ)

ggplot(data=obs.occ, aes(x=Dist_5km, y=ooC))+
   ggtitle("Observed occupancy of C given B")+
   geom_point()+
   geom_smooth(method="glm", method.args = list(family = quasibinomial(link = 'logit')))+
   facet_wrap(~ooC, ncol=1)

# and look at when the species are reversed
ggplot(data=obs.occ, aes(x=Dist_5km, y=ooB))+
   ggtitle("Observed occupancy of B given C")+
   geom_point()+
   geom_smooth(method="glm", method.args = list(family = quasibinomial(link = 'logit')))+
   facet_wrap(~ooB, ncol=1)






# cleanup
cleanup(ask=FALSE)

