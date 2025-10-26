# multi season multi-species model

# Fidino M, Simonis JL, Magle SB (2018) 
# A multi-species dynamic occupancy model to estimate local colonization-extinction 
# rates and patterns of co-occurrence between two or more interacting species. 
# Methods in Ecology and Evolution, online in advance of print. 
# https://doi.org/10.1111/2041-210x.13117
# 
# Several years of data on occurrences of 3 species. See their paper.
# 
# Data extracted from Dryad Digital Repository. 
# https://doi.org/10.5061/dryad.k5mp137

# 2018-12-16 Code contributed by Carl James Schwarz (cschwarz.stat.sfu.cs@gmail.com)


#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------

library(car)
library(ggplot2)
library(plyr)
library(readxl)
library(reshape2)
library(RMark)

# get the data
# importing detection / non-detection data for each of the 3 species at each camera trap.
# The data is the long format and we need to convert to the wide format

input.data <- read.csv(file.path("..","Fidino_MEE","data","fidino_sp_data.csv"),
                          as.is=TRUE, strip.white=TRUE, header=TRUE)
head(input.data)
input.data$Date <- as.Date(input.data$Date, "%m/%d/%Y")
input.data$Year <- as.numeric(format(input.data$Date, "%Y"))

# code the seasons numerically so sort properly. Notice the ordering
input.data$SeasonN <- car::recode(substr(input.data$Season,1,2),
                                  " 'SP'=2; 'SU'=3; 'FA'=4; 'WI'=1   ")
xtabs(~Date+SeasonWeek, data=input.data)
xtabs(~Year, data=input.data)

plyr::ddply(input.data, c("Year","SeasonN"), plyr::summarize,
            min.date=min(Date),
            max.date=max(Date),
            dif.date=max.date-min.date+1,
            n.date=length(unique(Date)))
head(input.data)

# According to the paper, they collapse each weeks of data to a single "visit" within a year-season,
# i.e. collapse the 28 days in a year-season to 4 visit in the year-season. Any detection on a camera
# during that period was declared a detection
input.data.red <- plyr::ddply(input.data, c("StationID","Year","SeasonN","Week"), plyr::summarize,
                              Coyote = any(Coyote),
                              Opossum= any(Opossum),
                              Raccoon= any(Raccoon))
head(input.data.red)


# How many Seasons are collected at each site.
any.data.seasons.station <- plyr::ddply(input.data.red, c("StationID","Year","SeasonN"), plyr::summarize, 
                                 any.data = sum(!is.na(Coyote)))
n.seasons.with.data      <- plyr::ddply(any.data.seasons.station, c("StationID"), plyr::summarize,
                                 n.seasons = sum(any.data>0))
xtabs(~n.seasons, n.seasons.with.data)

# According to the papery they exclude any station with less than 2 seasons of data reducing the data
# set to 103 records
exclude.stations <- n.seasons.with.data[ n.seasons.with.data$n.seasons < 2, c("StationID","n.seasons")]
exclude.stations

input.data.red <- input.data.red[ !input.data.red$StationID %in% exclude.stations$StationID,]

length(unique(input.data.red$Station))  # number of stations
nrow(unique(input.data.red[, c("Year","SeasonN","Week")]))     # = 13 (year-Seasons) x 4 weeks/season


coyote  <- reshape2::dcast(input.data.red, StationID ~ Year + SeasonN + Week, value.var="Coyote")
opossum <- reshape2::dcast(input.data.red, StationID ~ Year + SeasonN + Week, value.var="Opossum")
raccoon <- reshape2::dcast(input.data.red, StationID ~ Year + SeasonN + Week, value.var="Raccoon")

# MARK wants a 2 "digit" code xy for each visit where 
#     x = 0/1 if species A is not detected/detected
# and y = 0/1 if species B is not detected/detected

input.data2 <- 1*coyote[,-1] + 2*raccoon[,-1]
input.data2[1:2,]

input.chistory <- data.frame(lapply(as.data.frame(input.data2), as.character), stringsAsFactors=FALSE)

input.chistory <- data.frame(lapply(input.chistory, 
                           car::recode, " '0'='00'; '1'='10';'2'='01';'3'='11';", as.numeric=FALSE, as.factor=FALSE),
                           stringsAsFactors=FALSE)
input.history <- data.frame(ch=apply(input.chistory, 1, paste, sep="",collapse=""), freq=1)


# Change any NA to . in the chapter history
select <- grepl("NA", input.history$ch)
input.history[ select,][1:5,]

input.history$ch <- gsub("NA","..", input.history$ch, fixed=TRUE)
input.history[ select,][1:5,]




# Any site covariates?
site.covar <- read.csv(file.path("..","Fidino_MEE","data","fidino_covariate_data.csv"), header=TRUE, as.is=TRUE, strip.white=TRUE)
names(site.covar)
names(site.covar) <- make.names(names(site.covar))
names(site.covar)

site.covar <- site.covar[ order(site.covar$StationID),]

# standardized some variables
site.covar$House_std = scale(site.covar$House)

# Compute the first principal component on the site covariates (see paper)
URB <- princomp(site.covar[,c("House","Tree","Imp")], cor=TRUE, scores=TRUE)
URB$sdev^2 / sum(URB$sdev^2)

site.covar$UR <- URB$scores[,1]

ggplot(data=site.covar, aes(x=House, y=UR))+
   ggtitle("Relationship between Housing and UR index")+
   geom_point()

ggplot(data=site.covar, aes(x=Tree, y=UR))+
  ggtitle("Relationship between Tree and UR index")+
  geom_point()

head(site.covar)
# make sure that data is sorted properly and complete
setdiff(site.covar$StationID, coyote$StationID)
setdiff(coyote$StationID, site.covar$StationID)


# add the site covariate to the chapture hisotry
input.history <- cbind(input.history, site.covar)
head(input.history)

nrow(input.history)           # sites
nchar(input.history$ch[1])    # total number of visits 13 seasons x 4 visit/seasons x 2 digits/history element

time.intervals <-c( rep( c(rep(0, length(unique(input.data$Week))-1),1),  nrow(unique(input.data.red[, c("Year","SeasonN")]))-1 ),
                    rep(0, length(unique(input.data$Week))-1) )
length(time.intervals)



pcr.data <- process.data(data=input.history, time.intervals=time.intervals,
                         model="RD2SpGEConOcc")  # this parameterization is more stable

summary(pcr.data)

# are there any time varying covariates that need to be added to the ddls? 
pcr.ddl <- RMark::make.design.data(pcr.data)  # need for covariate predictions

# What are the parameter names for Single Season Single Species models
setup.parameters("RD2SpGEConOcc", check=TRUE)


#-----
# Use the psiAB parameterization 
#    Parameters of the model are 
#        PsiA - Initial probability that site is occupied by A only
#        PsiBA- Initial probability that site is occupied by B if A is present
#        PsiBa- Initial probability that site is occupied by B is A is absent
#
#    If species are independent thatn psiBA = psiBa.
#
#        GammaA  - Colonization by A given B absent at t-1
#        GammaB  - Colonization by B given A absent at t-1
#        GammaAB - Colonization by A given B present at t-1
#        GammaBA - Colonization by B given A present at t-1
 
#        EpsilonA  - extinction of A given B absent at t-1
#        EpsilonB  - extinction of B given A absent at t-1
#        EpsilonAB - extinction of A given B present at t-1
#        EpsilonBA - extinction of B given A present at t-1
#
#    This give the transition matrx (see (2) in the paper but notice that the 
#    From/To sides of the matrix have been transposed below
#                           To State
#                 U                            A                    B                          AB
#    From   U  (1-GammaA)(1-GammaB)      GammaA(1-GammaB)    (1-GammaA)GammaB             GammaA GammaB
#    State  A   EpsilonA (1-GammaBA) (1-EpsilonA)(1-GammaBA)   EpsilonGammaBA         (1-EpsilonA)(GammaBA)
#           B  (1-gammaAB)EpsilonB   GammaAB (EpsilonB)      (1-GammaAB)(1-EpsilonB)   GammaAB(1-EpsilonB)
#          AB   EpsilonAB EpsilonBA (1-EpsilonAB)EpsilonBA    EpsilonAB(1-EpsionBA)     (1-EpsionAB)(1-EpsilonBA)
#
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




# Get the list of models. NOtice NO equal signs here
model.list.csv <- textConnection("
PsiA,  PsiBA, PsiBa,   GammaA, GammaB, GammaAB, GammaBA,   EpsilonA, EpsilonB, EpsilonAB, EpsilonBA,    pA, pB, rA, rBA,rBa
~1,      ~1,    ~1,     ~1,      ~1,     ~1,     ~1,         ~1,       ~1,      ~1,         ~1,         ~1, ~1, ~1, ~1, ~1
~1,      ~1,    ~1,     ~1,      ~1,     ~1,     ~1,         ~1,       ~1,      ~1,         ~1,         ~1, ~1, ~1, ~1, ~SHARE
~UR,     ~UR,   ~UR,   ~UR,      ~UR,    ~UR,   ~UR,        ~UR,       ~UR,    ~UR,         ~UR,        ~1, ~1, ~1, ~1, ~1")



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
  
  # Now we need to set up the model parameters accounting for sharing
  model.parameters= list(
    PsiA     =list(formula=as.formula(eval(x$PsiA ))),
    PsiBA    =list(formula=as.formula(eval(x$PsiBA))),
    PsiBa    =list(formula=as.formula(eval(x$PsiBa))),
    GammaA   =list(formula=as.formula(eval(x$GammaA))),
    GammaB   =list(formula=as.formula(eval(x$GammaB))),
    GammaAB  =list(formula=as.formula(eval(x$GammaAB))),
    GammaBA  =list(formula=as.formula(eval(x$GammaBA))),
    EpsilonA =list(formula=as.formula(eval(x$EpsilonA))),
    EpsilonB =list(formula=as.formula(eval(x$EpsilonB))),
    EpsilonAB=list(formula=as.formula(eval(x$EpsilonAB))),
    EpsilonBA=list(formula=as.formula(eval(x$EpsilonBA))),
    pA       =list(formula=as.formula(eval(x$pA ))),
    pB       =list(formula=as.formula(eval(x$pB ))),
    rA       =list(formula=as.formula(eval(x$rA ))),
    rBA      =list(formula=as.formula(eval(x$rBA))),
    rBa      =list(formula=as.formula(eval(x$rBa)))
  )
  if(x$PsiBa == "~SHARE"){
    model.parameters$PsiBA  =list(formula=as.formula(eval(x$PsiBA)), share=TRUE)
    model.parameters$PsiBa  =NULL
  }
  if(x$rBa == '~SHARE'){
    model.parameters$rBA = list(formula=as.formula(eval(x$rBA)), share=TRUE)
    model.parameters$rBa = NULL
  }
  print(model.parameters)
  fit <- RMark::mark(input.data, ddl=input.ddl,
                     model="RD2SpGEConOcc",
                     model.parameters=model.parameters
                     #,brief=TRUE,output=FALSE, delete=TRUE
                     #,invisible=TRUE,output=TRUE  # set for debugging
  )
  mnumber <- paste("m...",formatC(x$model.number, width = 3, format = "d", flag = "0"),sep="")
  assign( mnumber, fit, envir=.GlobalEnv)
  #browser()
  fit
  
},input.data=pcr.data, input.ddl=pcr.ddl)


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

model.set <- RMark::collect.models( type="RD2SpGEConOcc")
model.set

names(model.set)
model.set$model.table



# make predictions as a function of urbanization
# We need to know the all.diff.index indices to the parameters
# Because we don't have any time effects in the model, we only need to find predictions for one session
setup.parameters("RD2SpGEConOcc", check=TRUE)


get.real(model.fits[[model.number]], "GammaA", se=TRUE)
GammaA <- RMark::covariate.predictions(model.set, indices=4, data=input.history)$estimates
GammaA$parameter <- 'GammaA'

get.real(model.fits[[model.number]], "GammaB", se=TRUE)
GammaB <- RMark::covariate.predictions(model.set, indices=16, data=input.history)$estimates
GammaB$parameter <- 'GammaB'

get.real(model.fits[[model.number]], "GammaAB", se=TRUE)
GammaAB <- RMark::covariate.predictions(model.set, indices=28, data=input.history)$estimates
GammaAB$parameter <- 'GammaAB'

get.real(model.fits[[model.number]], "GammaBA", se=TRUE)
GammaBA <- RMark::covariate.predictions(model.set, indices=40, data=input.history)$estimates
GammaBA$parameter <- 'GammaBA'


plotdata <- plyr::rbind.fill(GammaA, GammaB, GammaAB, GammaBA)
ggplot(data=plotdata, aes(x=UR, y=estimate))+
   ggtitle("Effects of urbanization on local colonization")+
   geom_point()+
   geom_ribbon(aes(ymin=lcl, ymax=ucl), alpha=0.2)+
   facet_wrap(~parameter, ncol=2)



# look at the extinction probabilities
get.real(model.fits[[model.number]], "EpsilonA", se=TRUE)
EpsilonA <- RMark::covariate.predictions(model.set, indices=52, data=input.history)$estimates
EpsilonA$parameter <- 'EpsilonA'

get.real(model.fits[[model.number]], "EpsilonB", se=TRUE)
EpsilonB <- RMark::covariate.predictions(model.set, indices=64, data=input.history)$estimates
EpsilonB$parameter <- 'EpsilonB'

get.real(model.fits[[model.number]], "EpsilonAB", se=TRUE)
EpsilonAB <- RMark::covariate.predictions(model.set, indices=76, data=input.history)$estimates
EpsilonAB$parameter <- 'EpsilonAB'

get.real(model.fits[[model.number]], "EpsilonBA", se=TRUE)
EpsilonBA <- RMark::covariate.predictions(model.set, indices=88, data=input.history)$estimates
EpsilonBA$parameter <- 'EpsilonBA'


plotdata <- plyr::rbind.fill(EpsilonA, EpsilonB, EpsilonAB, EpsilonBA)

ggplot(data=plotdata, aes(x=UR, y=estimate))+
  ggtitle("Effects of urbanization on local extinction")+
  geom_point()+
  geom_ribbon(aes(ymin=lcl, ymax=ucl), alpha=0.2)+
  facet_wrap(~parameter, ncol=2)




# cleanup
cleanup(ask=FALSE)

