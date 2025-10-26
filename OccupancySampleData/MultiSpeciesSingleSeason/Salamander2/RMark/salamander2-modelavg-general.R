# Multi species single season with a list of models

# Using RMark

# Co-occurence of Jordan's salamander (Plethodon jordani) (PJ) and
# members of Plethodon glutinosus (PG) in Great Smokey
# Mountains National Park (MacKenzie et al. 2004).

# 2018-09-27 Code contributed by Carl James Schwarz (cschwarz.stat.sfu.cs@gmail.com)


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
PG.data <- readxl::read_excel(file.path("..","Salamander2.xls"),
                              sheet="RawData", na='-',
                              col_names=FALSE,
                              range = "B3:F90")

PJ.data <- readxl::read_excel(file.path("..","Salamander2.xls"),
                              sheet="RawData", na='-',
                              col_names=FALSE,
                              range = "H3:L90")

# MARK wants a 2 "digit" code xy for each visit where 
#     x = 0/1 if species A is not detected/detected
# and y = 0/1 if species B is not detected/detected

input.data <- PG.data + 2*PJ.data
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
site.covar <- readxl::read_excel(file.path("..","Salamander2.xls"),
                                 sheet="RawData", col_names=TRUE,
                                 range = "N2:O90")
names(site.covar)
names(site.covar) <- make.names(names(site.covar))
names(site.covar)

# Standardize Elevation covariates
elevation.mean <- mean(site.covar$Elevation..m.)
elevation.std  <- sd  (site.covar$Elevation..m.)
site.covar$El <- (site.covar$Elevation..m. - elevation.mean)/elevation.std

# add the site covariate to the chapture hisotry
input.history <- cbind(input.history, site.covar)
head(input.history)

Nsites <- nrow(input.history)
Nvisits<- ncol(input.history)


sal.data <- process.data(data=input.history,
                         model="2SpecConOccup")  # this parameterization is more stable

summary(sal.data)

# If there are time-specific covariates, add them to the ddl's
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
PsiA,     PsiBA,   PsiBa,    pA,   pB,   rA,   rBA,   rBa
~1,        ~1,      ~1,     ~1,     ~1,  ~1,   ~1,    ~1
~1,        ~1,      ~1,     ~1,     ~1,  ~1,   ~1,    ~SHARE
~1,        ~1,      ~SHARE, ~1,     ~1,  ~1,   ~1,    ~1
~1,        ~1,      ~SHARE, ~1,     ~1,  ~1,   ~1,    ~SHARE
~El,       ~El,     ~El,    ~1,     ~1,  ~1,   ~1,    ~SHARE")


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
  
},input.data=sal.data,input.ddl=sal.ddl)


# examine individual model results
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

# make predictions as a function of elevation
# We need to know the indices to the parameters
# Look at the all.diff.index columns to get the index numbers

get.real(model.fits[[model.number]], "PsiA", se=TRUE)
PsiA.ma.e <- RMark::covariate.predictions(model.set, indices=1, data=input.history)$estimates
PsiA.ma.e$parameter <- 'PsiA'

get.real(model.fits[[model.number]], "PsiBa", se=TRUE)
PsiBa.ma.e <- RMark::covariate.predictions(model.set, indices=3, data=input.history)$estimates
PsiBa.ma.e$parameter <- 'PsiBa'

get.real(model.fits[[model.number]], "PsiBA", se=TRUE)
PsiBA.ma.e <- RMark::covariate.predictions(model.set, indices=2, data=input.history)$estimates
PsiBA.ma.e$parameter <- 'PsiBA'

# Not easy to get the se of dervied parameters at the covariate values
PsiB.ma.e <- data.frame(estimate=PsiBA.ma.e$estimate*PsiA.ma.e$estimate +
                                 PsiBa.ma.e$estimate*(1-PsiA.ma.e$estimate),
                        Elevation..m.=PsiA.ma.e$Elevation..m., stringAsFactors=FALSE)
PsiB.ma.e$parameter <- "PsiB"


plotdata <- plyr::rbind.fill(PsiA.ma.e,
                             PsiBA.ma.e,
                             PsiBa.ma.e,
                             PsiB.ma.e)

# It is clear that the fit is not satisfactory with a huge se for PsiBa and very steep decline.
# Try reversing the species and you get a much better fit
ggplot(data=plotdata, aes(x=Elevation..m., y=estimate))+
   ggtitle("Effects of elevation")+
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
obs.occ <- data.frame(Elevation..m. = site.covar$Elevation..m.,
                      ooPG = apply(PG.data,1,max, na.rm=TRUE),
                      ooPJ = apply(PJ.data,1,max, na.rm=TRUE))
head(obs.occ)

ggplot(data=obs.occ, aes(x=Elevation..m., y=ooPJ))+
   ggtitle("Observed occupancy of PJ given PG")+
   geom_point()+
   geom_smooth(method="glm", method.args = list(family = quasibinomial(link = 'logit')))+
   facet_wrap(~ooPG, ncol=1)

# and look at when the species are reversed
ggplot(data=obs.occ, aes(x=Elevation..m., y=ooPG))+
   ggtitle("Observed occupancy of PG given PJ")+
   geom_point()+
   geom_smooth(method="glm", method.args = list(family = quasibinomial(link = 'logit')))+
   facet_wrap(~ooPJ, ncol=1)






# cleanup
cleanup(ask=FALSE)

