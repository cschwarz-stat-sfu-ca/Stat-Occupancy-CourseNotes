# Brook Trout

# This is from MARK demo files on occupancy.
# Collected via electrofishing three 50 m sections of streams at 77 sites 
# in the Upper Chattachochee River basin. 
#       77 streams 3 occasions, 4 covariates: elevation, cross sectional area each occasion.

# Single Species Single Season Occupancy

# Fitting a single model using RPresence

# 2018-12-02 Code contributed by Neil Faught

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

library(car)
library(readxl)
library(RPresence)
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
input.history <- input.data # the history extracted
head(input.history)

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.history)
ncol(input.history)
range(input.history, na.rm=TRUE) # check that all values are either 0 or 1
sum(is.na(input.history))    # are there any missing values?


# Get the elevation information
elevation.data <- readxl::read_excel("../BrookTrout.xls",
                                     sheet="Elevation",
                                     na="-",
                                     col_names=TRUE)  
elevation.data$Elevation <- elevation.data$Elevation/1000  # standardized it a bit
elevation.data$Site      <- 1:nrow(elevation.data)
head(elevation.data)

# Get the cross sectional width
cross.data <- readxl::read_excel("../BrookTrout.xls",
                                 sheet="CrossSectionWidth",
                                 na="-",
                                 col_names=TRUE)  
head(cross.data)

# Convert to a survey covariate. You need to stack the data by columns
survey.cov <- data.frame(Site=1:nrow(cross.data),
                         visit=as.factor(rep(1:ncol(cross.data), each=nrow(cross.data))),
                         Cross =unlist(cross.data), stringsAsFactors=FALSE)

head(survey.cov)


# Create the *.pao file
trout.pao <- RPresence::createPao(input.history,
                                  unitcov=elevation.data,
                                  survcov=survey.cov,
                                  title='Brook Trout SSSS')
trout.pao

# Model where occupancy varies with elevation of stream
mod.elev <- RPresence::occMod(model=list(psi~Elevation, p~1), 
                            type="so", data=trout.pao)
summary(mod.elev)
head(mod.elev$real$psi)
mod.elev$real$p[seq(1, by=nrow(input.history), length.out=ncol(input.history)),]

# plot occupancy as a function of elevation of stream
plotdata = data.frame(mod.elev$real$psi, Elevation = elevation.data$Elevation)

ggplot(data=plotdata, aes(x=Elevation, y=est))+
  ggtitle("Occupancy as a function of elevation")+
  geom_point()+
  geom_ribbon(aes(ymin=lower_0.95, ymax=upper_0.95), alpha=0.2)+
  ylim(0,1)+
  ylab("Estimated occupancy")

