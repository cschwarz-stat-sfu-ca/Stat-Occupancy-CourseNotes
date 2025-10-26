# Brook Trout

# This is from MARK demo files on occupancy.
# Collected via electrofishing three 50 m sections of streams at 77 sites 
# in the Upper Chattachochee River basin. 
#       77 streams 3 occasions, 4 covariates: elevation, cross sectional area each occasion.

# Single Species Single Season Occupancy

# Fitting several models using RPresence

# 2018-08-15 Code contributed by Carl James Schwarz (cschwarz.stat.sfu.cs@gmail.com)

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



# define the list of models to fit
# Notice the commas between the column and the placement of the quotes
model.list.csv <- textConnection("
p,               psi
~1,              ~1
~1,              ~Elevation
~Cross,          ~1
~Cross,         ~Elevation")

model.list <- read.csv(model.list.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)
model.list


# fit the models
model.fits <- plyr::alply(model.list, 1, function(x,detect.pao){
  cat("\n\n***** Starting ", unlist(x), "\n")
  fit <- RPresence::occMod(model=list(as.formula(paste("psi",x$psi)),
                                      as.formula(paste("p"  ,x$p  ))),
                           data=detect.pao,type="so")
  fit
},detect.pao=trout.pao)




# Look the output from a specific model
check.model <- 1

names(model.fits[[check.model]])
model.fits[[check.model]]$beta

names(model.fits[[check.model]]$real)
model.fits[[check.model]]$real$psi[1:5,]
model.fits[[check.model]]$real$p[1:5,]

names(model.fits[[check.model]]$derived)
model.fits[[check.model]]$derived$psi_c[1:10,]
tail(model.fits[[check.model]]$derived$psi_c)




# Model averaging
aic.table <- RPresence::createAicTable(model.fits)
aic.table$table

names(aic.table)





# plot occupancy as a function of elevation of stream
psi.ma <- RPresence::modAvg(aic.table, param="psi")
head(psi.ma)
psi.ma$Site <- as.numeric(substring(row.names(psi.ma), 4+regexpr("unit",row.names(psi.ma), fixed=TRUE)))
plotdata <- merge(psi.ma, elevation.data)
head(plotdata)
ggplot(data=plotdata, aes(x=Elevation, y=est))+
   ggtitle("Occupancy as a function of elevation")+
   geom_point()+
   geom_ribbon(aes(ymin=lower_0.95, ymax=upper_0.95), alpha=0.2)+
   ylim(0,1)+
   ylab("Estimated occupancy")

# plot detection as a function of cross section
p.ma <- RPresence::modAvg(aic.table, param="p")
p.ma[grepl("unit1$", rownames(p.ma)),]
head(p.ma) 

p.ma$Site <- as.numeric(substring(row.names(p.ma), 4+regexpr("unit", row.names(p.ma), fixed=TRUE)))
p.ma$visit<- as.character(substr(row.names(p.ma), 2, -1+regexpr("_", row.names(p.ma), fixed=TRUE))) 

survey.cov$visit <- as.character(survey.cov$visit)

plotdata <- merge(p.ma, survey.cov)
head(plotdata)
ggplot(data=plotdata, aes(x=Cross, y=est))+
  ggtitle("Detection probability as a function of cross sectional width")+
  geom_point()+
  geom_ribbon(aes(ymin=lower_0.95, ymax=upper_0.95), alpha=.2)+
   ylim(0,1)+
   ylab("Estimated occupancy")

