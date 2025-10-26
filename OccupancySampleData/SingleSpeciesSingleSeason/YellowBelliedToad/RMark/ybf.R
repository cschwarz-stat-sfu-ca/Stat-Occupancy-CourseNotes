# Single Species Single Season Occupancy 

# Yellow-bellied toad
# Single Season Single Season occupancy 

#  RMark package

# 2018-11-15 Code contributed by Carl James Schwarz (cschwarz.stat.sfu.cs@gmail.com)

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

input.detect <- readxl::read_excel(file.path("..","YellowBelliedToad.xlsx"),
                                   sheet="detections")

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.detect)
ncol(input.detect)
range(input.detect, na.rm=TRUE)
sum(is.na(input.detect))

head(input.detect)


# Extract the history records and create a capture history
input.history <- data.frame(freq=1,
                            ch=apply(input.detect[,1:2],1,paste, collapse=""), stringsAsFactors=FALSE)
head(input.history)

# Change any NA to . in the chapter history
input.history$ch <- gsub("NA",".", input.history$ch, fixed=TRUE)
head(input.history)


# Get the site level covariates - none.



# Get the site x visit (sampling) covariates. 
jdate <- readxl::read_excel(file.path("..","YellowBelliedToad.xlsx"),
                            sheet="JulianDate")
range(jdate)

# standarize jdate
jdateS <- (jdate - 150)/10

names(jdate ) <- c('jdate1'   ,'jdate2')
names(jdateS) <- c('jdateS1'  ,'jdateS2')

# Sampling covariates must be added to the input.history as multiple columns
# We also need to create the quadratic terms in advance - this is a pain
jdateSQ <- jdateS
names(jdateSQ)<- c('jdateSQ1' , 'jdateSQ2')
jdateSQ <- jdateSQ^2


input.history <- cbind(input.history, jdate, jdateS, jdateSQ)
head(input.history)

ybf.data <- process.data(data=input.history,
                         model="Occupancy")
summary(ybf.data)

# Set up survey covariates in the ddl (in this case none)
ybf.ddl <- make.design.data(ybf.data)
ybf.ddl

# What are the parameter names for Single Season Single Species models
setup.parameters("Occupancy", check=TRUE)




####################################################################
# Fit a model
# Note that formula HAVE AN = SIGN
mod.fit <-  RMark::mark(ybf.data, ddl=ybf.ddl, 
                        model="Occupancy",
                        model.parameters=list(
                          Psi   =list(formula=~1),
                          p     =list(formula=~1)
                        )
)
summary(mod.fit)

# Look the objects returned in more details
names(mod.fit)
names(mod.fit$results)

# look at estimates on beta and original scale
mod.fit$results$beta  # on the logit scale

mod.fit$results$real# on the regular 0-1 scale for each site

# derived variabldes is the occupancy probability 
names(mod.fit$results$derived)

mod.fit$results$derived$Occupancy


# alternatively
get.real(mod.fit, "Psi", se=TRUE)
get.real(mod.fit, "Psi", pim=TRUE)


get.real(mod.fit, "p", se=TRUE)



#######################################
# Fit a model with p varying over time
# Note the use of the keyword "time" from the ddl which is a factor

ybf.ddl$p
str(ybf.ddl$p)

mod.fit2 <-  RMark::mark(ybf.data,  ddl=ybf.ddl,
                         model="Occupancy",
                         model.parameters=list(
                           Psi   =list(formula=~1),
                           p     =list(formula=~time)
                         )
)
summary(mod.fit2)

get.real(mod.fit2, "p", se=TRUE)

####################################################
# fit a model with p as a function of julian date (linear)
#

mod.fit3 <-  RMark::mark(ybf.data,  ddl=ybf.ddl,
                         model="Occupancy",
                         model.parameters=list(
                           Psi   =list(formula=~1),
                           p     =list(formula=~jdateS)
                         )
)
summary(mod.fit3)

# estimate of p at the "average" julian dates
get.real(mod.fit3, "p", se=TRUE) 

####################################################
# fit a model with p as a function of julian date (quadratinc)
#


mod.fit4 <-  RMark::mark(ybf.data,  ddl=ybf.ddl,
                         model="Occupancy",
                         model.parameters=list(
                           Psi   =list(formula=~1),
                           p     =list(formula=~jdateS+jdateSQ)
                         )
)
summary(mod.fit4)

# estimate of p at the "average" julian dates and average of the quadratic term - this is silly
get.real(mod.fit4, "p", se=TRUE) 

##############################################
# Collect models
model.set <- RMark::collect.models( type="Occupancy")
model.set


names(model.set)
model.set$model.table


# model averaged values
get.real(mod.fit , "Psi", se=TRUE)
get.real(mod.fit2, "Psi", se=TRUE)
get.real(mod.fit3, "Psi", se=TRUE)

Psi.ma <- RMark::model.average(model.set, param="Psi")
Psi.ma



# It is often convenient to estimate parameters (such as p) for each site using the
# values of the covariates specific for that site.
# This is similar to RPresence which gives estimates at each site.

# We create a data frame with the covariates etc as measured on the input.history dataframe.
# But we also need to add the appropriate index number for psi to the covariate data frame to account
# for the groups definition

p.covariates <- input.history
p.covariates$freq <- NULL
p.covariates$ch   <- NULL

ybf.ddl$p
p.covariates <- merge(ybf.ddl$p, p.covariates)
head(p.covariates)

# we only want a prediction for that particular model.index.
# we need to create a variable "index" that has the model.index
p.covariates$index <- p.covariates$model.index

p.pred <- covariate.predictions(model.set,
                                   data=p.covariates)

# Notice that there is problem in that if model.index =1, it used jdate1 as the predictor
# and if model.index==2, it used jdate2 as the predictor. There is no easy around thi!
head(p.pred$estimates)


ggplot(data=NULL, aes( y=estimate))+
  ggtitle("Model averaged predictions of detection for each site")+
  geom_point(data=p.pred$estimates[ p.pred$estimates$index==2,], aes(x=jdate2))+
  geom_point(data=p.pred$estimates[ p.pred$estimates$index==1,], aes(x=jdate1))+
  geom_ribbon(data=p.pred$estimates[ p.pred$estimates$index==1,],aes(x=jdate1, ymin=lcl, ymax=ucl), alpha=0.2)+
  geom_ribbon(data=p.pred$estimates[ p.pred$estimates$index==2,],aes(x=jdate2, ymin=lcl, ymax=ucl), alpha=0.2)


# cleanup
cleanup(ask=FALSE)
