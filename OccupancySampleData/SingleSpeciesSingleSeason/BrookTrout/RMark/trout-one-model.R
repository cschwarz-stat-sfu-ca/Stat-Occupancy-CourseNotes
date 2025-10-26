# Brook Trout

# This is from MARK demo files on occupancy.
# Collected via electrofishing three 50 m sections of streams at 77 sites 
# in the Upper Chattachochee River basin. 
#       77 streams 3 occasions, 4 covariates: elevation, cross sectional area each occasion.

# Single Species Single Season Occupancy


# Fitting models using RMark

# 2018-11-09 Code contributed by Neil Faught

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

input.data <- readxl::read_excel("../BrookTrout.xls",
                                 sheet="DetectionHistory",
                                 na="-",
                                 col_names=FALSE)  # notice no column names in row 1 of data file. 

head(input.data)


# Extract the history records

# Extract the history records and create a capture history
input.history <- data.frame(freq=1,
                            ch=apply(input.data,1,paste, collapse=""), stringsAsFactors=FALSE)
head(input.history)

# Change any NA to . in the chapter history
input.history$ch <- gsub("NA",".", input.history$ch, fixed=TRUE)
head(input.history)

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.history)




# Get the pond information
# Get the elevation information
elevation.data <- readxl::read_excel("../BrookTrout.xls",
                                     sheet="Elevation",
                                     na="-",
                                     col_names=TRUE)  
input.history$Elevation <- elevation.data$Elevation/1000  # standardized it a bit
head(input.history)




# Get the cross sectional width
cross.data <- readxl::read_excel("../BrookTrout.xls",
                                 sheet="CrossSectionWidth",
                                 na="-",
                                 col_names=TRUE)  
head(cross.data)
# Notice that RMark does not allow missing values in time-varying covariates, even when visits are not made
# so do not set these to missing as in RPresence

input.history <- cbind(input.history, cross.data)
head(input.history)

# Create the data structure
trout.data <- process.data(data=input.history,
                           model="Occupancy")
summary(trout.data)

# If survey covariates are present, modify the ddl.
trout.ddl <- make.design.data(trout.data)
trout.ddl

# What are the parameter names for Single Season Single Species models
setup.parameters("Occupancy", check=TRUE)

# Fit a model
# Models where detection probability (p) varies with temperature (Cross), and
# Occupancy (Psi) varies with Elevation
# Notice that RMark does not allow missing values in time-varying covariates, even when visits are not made
mod.fit <-  RMark::mark(trout.data, ddl = trout.ddl,
                        model="Occupancy",
                        model.parameters=list(
                          Psi   =list(formula=~Elevation),
                          p     =list(formula=~Cross) 
                        )
)

summary(mod.fit)


# Look the objects returned in more details
names(mod.fit)

# look at estimates on beta and original scale
mod.fit$results$beta  # on the logit scale

mod.fit$results$real# on the regular 0-1 scale for each site

# derived variables is the occupancy probability 
names(mod.fit$results$derived)

mod.fit$results$derived$Occupancy


# get the psi-values as function of the covariates
trout.ddl$Psi # see the index numbers

Elev.df <-data.frame(Elevation=seq(min(input.history$Elevation,na.rm=TRUE),
                                   max(input.history$Elevation, na.rm=TRUE),.1)) 
pred.psi <- covariate.predictions(mod.fit, indices=4, data=Elev.df)
head(pred.psi$estimates)

ggplot(data=pred.psi$estimates, aes(x=covdata, y=estimate))+
  ggtitle("Occupancy probability as a function of elevation")+
  geom_point()+
  geom_ribbon(aes(ymin=lcl, ymax=ucl), alpha=.2)+
  ylim(0,1)+
  ylab("Estimated probability of occupancy")+
  xlab("Elevation (km)")


# make a plot of the probability of detection as a function of Cross sectional width
Cross.df <- data.frame(Cross1=seq(min(cross.data,na.rm=TRUE),max(cross.data, na.rm=TRUE),0.1))

trout.ddl$p # see the index numbers

pred.p <- covariate.predictions(mod.fit, indices=1, data=Cross.df)
head(pred.p$estimates)

plotdata <- cbind(Cross.df, pred.p$estimates)
head(plotdata)

ggplot(data=plotdata, aes(x=Cross1, y=estimate))+
  ggtitle("Detection probability as a function of cross sectional width")+
  geom_point()+
  geom_ribbon(aes(ymin=lcl, ymax=ucl), alpha=.2)+
  ylim(0,1)

# It is often convenient to estimate occupancy for each site using the
# values of the covariates specific for that site.
# This is similar to RPresence which gives estimates at each site.

# We create a data frame with the covariates etc as measure on the input.history dataframe.
# But we also need to add the appropriate index number for psi to the covariate data frame to account
# for the groups definition

site.covariates <- input.history
site.covariates$freq <- NULL
site.covariates$ch   <- NULL

trout.ddl$Psi
site.covariates <- merge(trout.ddl$Psi, site.covariates)
head(site.covariates)

# we only want a prediction for that particular model.index.
# we need to create a variable "index" that has the model.index
site.covariates$index <- site.covariates$model.index

site.pred <- covariate.predictions(mod.fit,
                                   data=site.covariates)
ggplot(site.pred$estimates, aes(x=Elevation, y=estimate))+
  ggtitle("Model predictions of occupancy for each site")+
  geom_point()+
  geom_ribbon(aes(ymin=lcl, ymax=ucl), alpha=0.2)

# cleanup
cleanup(ask=FALSE)
