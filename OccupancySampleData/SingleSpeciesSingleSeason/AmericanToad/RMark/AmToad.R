# American Toad

# Extracted from
#    Darryl I. MacKenzie, et al. 2002.
#    Estimating site occupancy rates when detection probabilities are less than one.
#    Ecology 83:2248-2255.
#    doi:10.1890/0012-9658(2002)083[2248:ESORWD]2.0.CO;2]

# 29 sites with 82 sampling occasions in 2000.
# Volunteers visited sites and recorded presence/absence of toads by calls.
# Habitat (type of pond, permanent or ephemeral) and temperature at visit recorded.

# Single Species Single Season Occupancy

# Fitting models using RMark


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

input.data <- readxl::read_excel("../AmericanToad.xls",
                                 sheet="AmToadDetectionHistories",
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
pond.data <- readxl::read_excel("../AmericanToad.xls",
                                sheet="Pond",
                                na="-",
                                col_names=TRUE)  
head(pond.data)
input.history$Pond <- as.factor(car::recode(pond.data$Pond,
                                        "1='P'; 0='E'; " ))
head(input.history)




# Get the temperature data
temp.data <- readxl::read_excel("../AmericanToad.xls",
                                sheet="AmToadTemperature",
                                na="-",
                                col_names=FALSE)  # notice no column names in row 1 of data file. 
head(temp.data)
# Notice that RMark does not allow missing values in time-varying covariates, even when visits are not made
# so do not set these to missing as in RPresence
# temp.data[ is.na(input.data)] <- NA


colnames(temp.data) = paste("Temp", 1:ncol(temp.data), sep="") # assign the observer at each time point
head(temp.data)

input.history <- cbind(input.history, temp.data)
head(input.history)

# Illustration of using categorical covariate in the modelling process by creating groups
amtoad.data <- process.data(data=input.history, group="Pond",
                         model="Occupancy")
summary(amtoad.data)

# If survey covariates are present, modify the ddl
amtoad.ddl <- make.design.data(amtoad.data)
amtoad.ddl

# What are the parameter names for Single Season Single Species models
setup.parameters("Occupancy", check=TRUE)


# Fit a model
# Notice that RMark does not allow missing values in time-varying covariates, even when visits are not made

# Note that formula have an equal sign
mod.fit <-  RMark::mark(amtoad.data, ddl=amtoad.ddl,
                        model="Occupancy",
                        model.parameters=list(
                          Psi   =list(formula=~Pond),
                          p     =list(formula=~Temp) # need to specify rest of obsevers after the intercept 
                        )
)
summary(mod.fit)

# Look the objects returned in more details
names(mod.fit)
names(mod.fit$results)

# look at estimates on beta and original scale
mod.fit$results$beta  # on the logit scale

mod.fit$results$real# on the regular 0-1 scale for each site

# derived variables is the occupancy probability 
names(mod.fit$results$derived)

mod.fit$results$derived$Occupancy


# get the two psi values and their covariance
get.real(mod.fit, "Psi", se=TRUE, vcv=TRUE)

get.real(mod.fit, "Psi", pim=TRUE)


# make a plot of the probability of detection as a function of temperature
Temp.df <- data.frame(Temp1=seq(min(temp.data,na.rm=TRUE),max(temp.data, na.rm=TRUE),1))

amtoad.ddl$p # see the index numbers

pred.p <- covariate.predictions(mod.fit, indices=1, data=Temp.df)
head(pred.p$estimates)

plotdata <- cbind(Temp.df, pred.p$estimates)
head(plotdata)

ggplot(data=plotdata, aes(x=Temp1, y=estimate))+
      ggtitle("Detection probability as a function of temperature")+
      geom_point()+
      geom_ribbon(aes(ymin=lcl, ymax=ucl), alpha=.2)+
      ylim(0,1)

# covariate predictions for categorical covariates
# Occupancy predictions for different pond types
amtoad.ddl$Psi
fit.psi = tail(mod.fit$results$real,2)
fit.psi$pond = c("E","P")
fit.psi

ggplot(data=fit.psi, aes(x=pond, y=estimate))+
  geom_point()+
  geom_errorbar(aes(ymin=lcl, 
                    ymax=ucl, width=0.2))+
  ylim(0,1)

# cleanup
cleanup(ask=FALSE)


