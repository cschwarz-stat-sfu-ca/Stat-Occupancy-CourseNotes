# Single Species Single Season Occupancy models using RMark 

# Blue Gross Beaks.
# Downloaded from https://sites.google.com/site/asrworkshop/home/schedule/r-occupancy-1
# Changed BQI to Y/N from 0/1

# An occupancy study was made on Blue Grosbeaks (Guiraca caerulea) 
# on 41 old fields planted to longleaf pines (Pinus palustris) 
# in southern Georgia, USA. 

# Surveys were 500 m transects across each field 
# and were completed three times during the breeding season in 2001.

# Columns in the file are:
#    field - field number
#    v1, v2, v3 -  detection histories for each site on each of 3 visit during the 2001 breeding season.    
#    bqi - bobwhite quality initiative applied to this field
#    field.size - size of the files
#    crop.hist - crop history
#    crop1, crop2 - indicator variables for the crop history
#    count1, count2, count3 - are actual counts of birds detected in each visit

# 2018-08-15 Code contributed by Carl James Schwarz (cschwarz.stat.sfu.cs@gmail.com)



#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

#  RMark package

library(readxl)
library(RMark)
library(ggplot2)

# get the data read in
# Data for detections should be a data frame with rows corresponding to sites
# and columns to visits.
# The usual 1=detected; 0=not detected; NA=not visited is used.

input.data <- read.csv(file.path("..","blgr.csv"), 
                       header=TRUE, as.is=TRUE, strip.white=TRUE) 
head(input.data)

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.data)
range(input.data[, c("v1","v2","v3")], na.rm=TRUE)
sum(is.na(input.data[, c("v1","v2","v3")]))

input.history <- input.data[, c("v1","v2","v3")]
head(input.history)

site.covar <- input.data[, c("field","field.size","bqi")]
site.covar$logFS <- log(site.covar$field.size)
head(site.covar)

# Extract the history records and create a capture history
input.history <- data.frame(freq=1,
                            ch=apply(input.history,1,paste, collapse=""), 
                            stringsAsFactors=FALSE)
head(input.history)

# Change any NA to . in the capture history
input.history$ch <- gsub("NA",".", input.history$ch, fixed=TRUE)
head(input.history)

#Add field size data to the capture history
input.history = cbind(input.history, site.covar)
head(input.history)

#Create the data structure
grossbeak.data <- process.data(data = input.history,
                               group="bqi",
                               model = "Occupancy")
summary(grossbeak.data)

# modify the ddl if needed (e.g. for site level covariates)
grossbeak.ddl <- make.design.data(grossbeak.data)
grossbeak.ddl

# What are the parameter names for Single Season Single Species models
setup.parameters("Occupancy", check=TRUE)




# Get the list of models. NOtice NO equal signs here
model.list.csv <- textConnection("
                                 p,            Psi
                                 ~1,         ~1
                                 ~time,      ~1
                                 ~1,         ~logFS
                                 ~1,         ~bqi")

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
  
},input.data=grossbeak.data, input.ddl=grossbeak.ddl)


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





# covariate predictions for continuous covariates
# such as logFS

# need to find the index number of the parameter
grossbeak.ddl$Psi

new.logFS <- data.frame(logFS=seq(min(site.covar$logFS),
                                  max(site.covar$logFS), .05))
psi.ma <- covariate.predictions(model.set,
                                indices=7,
                                data=new.logFS)
psi.ma$estimates

ggplot(data=psi.ma$estimates, aes(x=covdata, y=estimate))+
   geom_point()+
   geom_ribbon(aes(ymin=lcl, ymax=ucl), alpha=0.2)



# covariate predictions for categorical covariates
# such a bqi
grossbeak.ddl$Psi

psi.ma <- model.average(model.set, param="Psi", vcv=TRUE)
psi.ma$estimates


ggplot(data=psi.ma$estimates, aes(x=bqi, y=estimate))+
  geom_point()+
  geom_errorbar(aes(ymin=lcl, 
                    ymax=ucl), width=0.2)+
  ylim(0,1)


# or use covariate predictions
psi.ma2 <- covariate.predictions(model.set, indices=7:8)
psi.ma2$estimates <- cbind(psi.ma2$estimates, grossbeak.ddl$Psi[,"bqi", drop=FALSE])
psi.ma2$estimates

ggplot(data=psi.ma2$estimates, aes(x=bqi, y=estimate))+
  geom_point()+
  geom_errorbar(aes(ymin=lcl, 
                    ymax=ucl), width=0.2)+
  ylim(0,1)

# caution when you have both continous and categorical covariates
# as it is not clear what value of the continuous covariate is used
psi.ma
psi.ma2$estimates

# It is often convenient to estimate occupancy for each site using the
# values of the covariates specific for that site.
# This is similar to RPresence which gives estimates at each site.

# We create a data frame with the covariates etc as measured on the input.history dataframe.
# But we also need to add the appropriate index number for psi to the covariate data frame to account
# for the groups definition

site.covariates <- input.history
site.covariates$freq <- NULL
site.covariates$ch   <- NULL

grossbeak.ddl$Psi
site.covariates <- merge(grossbeak.ddl$Psi, site.covariates)
head(site.covariates)

# we only want a prediction for that particular model.index.
# we need to create a variable "index" that has the model.index
site.covariates$index <- site.covariates$model.index

site.pred <- covariate.predictions(model.set,
                                   data=site.covariates)

ggplot(site.pred$estimates, aes(x=logFS, y=estimate))+
  ggtitle("Model averaged predictions of occupancy for each site")+
  geom_point()+
  geom_ribbon(aes(ymin=lcl, ymax=ucl), alpha=0.2)+
  facet_wrap(~bqi, ncol=1)

#cleanup
cleanup(ask=FALSE)

