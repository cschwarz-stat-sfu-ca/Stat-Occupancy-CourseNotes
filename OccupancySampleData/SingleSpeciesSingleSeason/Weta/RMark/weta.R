# Single Species Single Season Occupancy 

# Weta Example

#  RMark -
#      categorical variable for site level covariate can be modelled directly as a covariate
#         or as a grouping variable with the same results.
#      showing how to use a covariate that varies by sites x visit (e.g. observers)


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

input.data <- readxl::read_excel(file.path("..","weta.xls"),
                                    sheet="detection_histories",
                                    na="-",
                                    col_names=FALSE)  # notice no column names in row 1 of data file. 

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.data)
ncol(input.data)
range(input.data, na.rm=TRUE)
sum(is.na(input.data))

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




# Get the site level covariates
site_covar <- readxl::read_excel(file.path("..","weta.xls"),
                                 sheet="site_covar",
                                 na="-",
                                 col_names=TRUE)  # notice col_names in row 1 of table. 
head(site_covar)

# Create a site level covariate that is a categorical variable rather 
# than indicator variables. Be sure that it is declared as a FACTOR
input.history$BrowCat <- factor(site_covar$BrowseCat)

head(input.history)


# Get the individual covariates. 
obs1 <- readxl::read_excel(file.path("..","weta.xls"),
                           sheet="Obs1",
                           na="-",
                           col_names=FALSE) 
obs2 <- readxl::read_excel(file.path("..","weta.xls"),
                           sheet="Obs2",
                           na="-",
                           col_names=FALSE) 
obs3 <- readxl::read_excel(file.path("..","weta.xls"),
                           sheet="Obs3",
                           na="-",
                           col_names=FALSE) 

Obs <-obs1*1 + obs2*2 + obs3*3
Obs<- as.data.frame(lapply(Obs, as.character), stringsAsFactors=FALSE)
Obs[ is.na(Obs)] <- "."  # set all missing values to .
head(Obs)

# make each column a factor variable. Be sure to use 
Obs[] <- lapply(Obs, as.factor)
str(Obs)

colnames(Obs) = paste("Observer", 1:5, sep="") # assign the observer at each time point
head(Obs)

Obs.df <- make.time.factor(Obs, var.name="Observer", times=1:5, intercept=1)
head(Obs.df) # Column names are indicator variables for Obsever x at time y

input.history <- cbind(input.history, Obs.df)
head(input.history)


# Illustration of using categorical covariate in the modelling process by creating groups
weta.data <- process.data(data=input.history, group="BrowCat",
                         model="Occupancy")
summary(weta.data)

# What are the parameter names for Single Season Single Species models
setup.parameters("Occupancy", check=TRUE)

# modify the ddl if needed (e.g. for visit level covariates) - none in this case
weta.ddl <- make.design.data(weta.data)
weta.ddl

#----------------------------------------------------------------
# some simple models
mod.fit1 <-  RMark::mark(weta.data, ddl=weta.ddl,
                        model="Occupancy",
                        model.parameters=list(
                          Psi   =list(formula=~1),
                          p     =list(formula=~1) 
                        )
)

summary(mod.fit1)

get.real(mod.fit1, "Psi", se=TRUE)

get.real(mod.fit1, "p",   se=TRUE)




#----------------------------------------------------------------
# some simple models
mod.fit2 <-  RMark::mark(weta.data, ddl = weta.ddl,
                         model="Occupancy",
                         model.parameters=list(
                           Psi   =list(formula=~BrowCat),
                           p     =list(formula=~1) 
                         )
)

summary(mod.fit2)

get.real(mod.fit2, "Psi", se=TRUE)

get.real(mod.fit2, "p",   se=TRUE)



#----------------------------------------------------------------
# some simple models
mod.fit3 <-  RMark::mark(weta.data, ddl=weta.ddl,
                         model="Occupancy",
                         model.parameters=list(
                           Psi   =list(formula=~BrowCat),
                           p     =list(formula=~BrowCat) 
                         )
)

summary(mod.fit3)

get.real(mod.fit3, "Psi", se=TRUE)

get.real(mod.fit3, "p",   se=TRUE)


#----------------------------------------------------------------
# some simple models

# Fit a model
# Note that formula have an equal sign
mod.fit4 <-  RMark::mark(weta.data, ddl=weta.ddl,
                        model="Occupancy",
                        model.parameters=list(
                          Psi   =list(formula=~BrowCat),
                          p     =list(formula=~Observer2+Observer3) # need to specify rest of obsevers after the intercept 
                        )
)
summary(mod.fit4)

# Look the objects returned in more details
names(mod.fit4)
names(mod.fit4$results)

# look at estimates on beta and original scale
mod.fit4$results$beta  # on the logit scale

mod.fit4$results$real# on the regular 0-1 scale for each site

# derived variabldes is the occupancy probability 
names(mod.fit4$results$derived)

mod.fit4$results$derived$Occupancy


# get the two psi values and their covariance
get.real(mod.fit4, "Psi", se=TRUE, vcv=TRUE)

get.real(mod.fit4, "Psi", pim=TRUE)




# Estimate of p are much ore complicated when you have site x visit covariates
# These don't make much sense as they are computed at the average
# value of the observer for time 1:5
get.real(mod.fit4, "p", se=TRUE)

# get the estimated directly from the beta matrix
mod.fit4$results$beta  # on the logit scale
obs.logit <- matrix(c(1,0,0,0,0,
                      1,1,0,0,0,
                      1,0,1,0,0), ncol=5, byrow=TRUE) %*% mod.fit4$results$beta$estimate
expit <- function(x){1/(1+exp(-x))}
expit(obs.logit)

# Or we can get covariate predictions
weta.ddl$p # see the index numbers
# use the model.index index nubers

# create a special matrix, or use the data matrix for each site x visit

obs.df <-data.frame(Observer21=c(0,1,0),
                    Observer22=c(0,1,0),
                    Observer23=c(0,1,0),
                    Observer24=c(0,1,0),
                    Observer25=c(0,1,0),
                    Observer31=c(0,0,1),
                    Observer32=c(0,0,1),
                    Observer33=c(0,0,1),
                    Observer34=c(0,0,1),
                    Observer35=c(0,0,1))
 
obs.p <- covariate.predictions(mod.fit4, indices=1:10, data=obs.df)
obs.p$estimates




##----------------------------------------------------------------------------
mod.fit5 <-  RMark::mark(weta.data, ddl = weta.ddl,
                         model="Occupancy",
                         model.parameters=list(
                           Psi   =list(formula=~BrowCat),
                           p     =list(formula=~time+Observer2+Observer3) # need to specify rest of obsevers after the intercept 
                         )
)
summary(mod.fit5)


##----------------------------------------------------------------------------
mod.fit6 <-  RMark::mark(weta.data, ddl = weta.ddl,
                         model="Occupancy",
                         model.parameters=list(
                           Psi   =list(formula=~BrowCat),
                           p     =list(formula=~time+BrowCat+Observer2+Observer3) # need to specify rest of obsevers after the intercept 
                         )
)
summary(mod.fit6)

##----------------------------------------------------------------------------
mod.fit7 <-  RMark::mark(weta.data, ddl = weta.ddl,
                         model="Occupancy",
                         model.parameters=list(
                           Psi   =list(formula=~1),
                           p     =list(formula=~time+BrowCat+Observer2+Observer3) # need to specify rest of obsevers after the intercept 
                         )
)
summary(mod.fit7)

##----------------------------------------------------------------------------
mod.fit8 <-  RMark::mark(weta.data, ddl = weta.ddl,
                         model="Occupancy",
                         model.parameters=list(
                           Psi   =list(formula=~1),
                           p     =list(formula=~time+Observer2+Observer3) # need to specify rest of obsevers after the intercept 
                         )
)
summary(mod.fit8)



######################################################333
######################################################333
######################################################333

model.set <- collect.models(type="Occupancy")
model.set


psi.ma <- RMark::model.average(model.set, parameter="Psi", vcv=TRUE)
names(psi.ma)
head(psi.ma$estimates)


ggplot(data=psi.ma$estimates, aes(x=BrowCat, y=estimate))+
  geom_point()+
  geom_errorbar(aes(ymin=lcl, 
                    ymax=ucl), width=0.2)+
  ylim(0,1)+
  ylab("Occupancy (with 95% ci)")+xlab("Browse Category")

# or use covariate predictions to estimate the occupancy at different
# browse categories

# use the get.real to figure out the index numbers
# by looking at the alt.diff.index numbers
get.real(mod.fit8, parameter="Psi", se=TRUE)

psi.ma2 <- covariate.predictions(model.set, indices=11:12)
head(psi.ma2$estimates)
psi.ma2$estimates <- merge(psi.ma2$estimates, 
                           weta.ddl$Psi[,c("model.index","BrowCat")],
                           by.x="par.index", by.y="model.index")
psi.ma2$estimates

ggplot(data=psi.ma2$estimates, aes(x=BrowCat, y=estimate))+
  geom_point()+
  geom_errorbar(aes(ymin=lcl, 
                    ymax=ucl), width=0.2)+
  ylim(0,1)+
  ylab("Occupancy (with 95% ci)")+xlab("Browse Category")


# Much more difficult to compute model averaged values of p

# cleanup
cleanup(ask=FALSE)


