# Example using model averaging on the pronghorn dataset
# 256 sites, two sampling occasions per site
# Covariates used in this example:
# 1. sagebrush (continuos) - Sagebrush density
# 2. aspect (Categorical) - Compass direction slope faces (N,S,E,W)

# Code contributed by Neil Faught - 26/11/2018

#  RMark package

library(readxl)
library(RMark)
library(ggplot2)

# get the data read in
# Data for detections should be a data frame with rows corresponding to sites
# and columns to visits.
# The usual 1=detected; 0=not detected; NA=not visited is used.

input.data<-read.csv(file.path("..","pronghorn.csv"), header=TRUE, as.is=TRUE, strip.white=TRUE)
head(input.data)

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.data)
range(input.data[, c("Survey.1","Survey.2")], na.rm=TRUE)
sum(is.na(input.data[, c("Survey.1","Survey.2")]))

input.history <- input.data[, c("Survey.1","Survey.2")]
head(input.history)

site.covar <- input.data[, c("sagebrush","aspect")]
head(site.covar)
#Convert aspect to a factor
site.covar$aspect = as.factor(site.covar$aspect)

# Extract the history records and create a capture history
input.history <- data.frame(freq=1,
                            ch=apply(input.history,1,paste, collapse=""), 
                            stringsAsFactors=FALSE)
head(input.history)

# Change any NA to . in the capture history
input.history$ch <- gsub("NA",".", input.history$ch, fixed=TRUE)
head(input.history)

#Add aspect and sagebrush data to the capture history
input.history = cbind(input.history, site.covar)
head(input.history)

#Create the data structure
pronghorn.data <- process.data(data = input.history,
                               group="aspect",
                               model = "Occupancy")
summary(pronghorn.data)

# modify the ddl if needed (e.g. for survey/visit level covariates)
pronghorn.ddl <- make.design.data(pronghorn.data)
pronghorn.ddl

# What are the parameter names for Single Season Single Species models
setup.parameters("Occupancy", check=TRUE)

# Get the list of models. NOtice NO equal signs here
model.list.csv <- textConnection("
                                 p,            Psi
                                 ~1,         ~1
                                 ~time,      ~1
                                 ~1,         ~aspect
                                 ~1,         ~sagebrush")

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
  
},input.data=pronghorn.data, input.ddl=pronghorn.ddl)

# examine individula model results
model.number <-4

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
get.real(model.set[[4]], "Psi", se=TRUE)


# covariate predictions for continuous covariates
# such as sagebrush

# need to find the index number of the parameter
pronghorn.ddl$Psi

new.sagebrush <- data.frame(sagebrush=seq(min(site.covar$sagebrush),
                                  max(site.covar$sagebrush), .05))
# Note that this just makes predictions for observations from the east aspect
# Change the indices= line to view plots for different aspects. To view plots
# from all aspects at once, use indices=9:12
psi.ma <- covariate.predictions(model.set,
                                indices=9,
                                data=new.sagebrush)
psi.ma$estimates

ggplot(data=psi.ma$estimates, aes(x=covdata, y=estimate))+
  geom_point()+
  geom_ribbon(aes(ymin=lcl, ymax=ucl), alpha=0.2)+
  labs(x = "Sagebrush", y = "Estimated Occupancy")

# covariate predictions for categorical covariates
# such as aspect
pronghorn.ddl$Psi

psi.ma <- model.average(model.set, param="Psi", vcv=TRUE)
names(psi.ma)
psi.ma$estimates

ggplot(data=psi.ma$estimates, aes(x=aspect, y=estimate))+
  geom_point()+
  geom_errorbar(aes(ymin=lcl, 
                    ymax=ucl), width=0.2)+
  ylim(0,1)+
  labs(x = "Aspect", y = "Estimated Occupancy")

# or use covariate predictions
psi.ma2 <- covariate.predictions(model.set, indices=9:12)
psi.ma2$estimates <- cbind(psi.ma2$estimates, pronghorn.ddl$Psi[,"aspect", drop=FALSE])
psi.ma2$estimates

ggplot(data=psi.ma2$estimates, aes(x=aspect, y=estimate))+
  geom_point()+
  geom_errorbar(aes(ymin=lcl, 
                    ymax=ucl), width=0.2)+
  ylim(0,1)+
  labs(x = "Aspect", y = "Estimated Occupancy")


# caution when you have both continous and categorical covariates
# as it is not clear what value of the continuous covariate is used
psi.ma
psi.ma2$estimates

# It is often convenient to estimate occupancy for each site using the
# values of the covariates specific for that site.
# This is similar to RPresence which gives estimates at each site.

# We create a data frame with the covariates etc as measure on the input.history dataframe.
# But we also need to add the appropriate index number for psi to the covariate data frame to account
# for the groups definition

site.covariates <- input.history
site.covariates$freq <- NULL
site.covariates$ch   <- NULL

pronghorn.ddl$Psi
site.covariates <- merge(pronghorn.ddl$Psi, site.covariates)
head(site.covariates)

# we only want a prediction for that particular model.index.
# we need to create a variable "index" that has the model.index
site.covariates$index <- site.covariates$model.index

site.pred <- covariate.predictions(model.set,
                                   data=site.covariates)

ggplot(site.pred$estimates, aes(x=sagebrush, y=estimate))+
  ggtitle("Model averaged predictions of occupancy for each site")+
  geom_point()+
  geom_ribbon(aes(ymin=lcl, ymax=ucl), alpha=0.2)+
  facet_wrap(~aspect, ncol=1)

#cleanup
cleanup(ask=FALSE)
