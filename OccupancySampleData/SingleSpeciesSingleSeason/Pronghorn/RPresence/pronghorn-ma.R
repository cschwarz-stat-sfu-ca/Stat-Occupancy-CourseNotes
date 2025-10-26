# Example using model averaging on the pronghorn dataset
# 256 sites, two sampling occasions per site
# Covariates used in this example:
# 1. sagebrush (continuos) - Sagebrush density
# 2. aspect (Categorical) - Compass direction slope faces (N,S,E,W)

# 2018-11-26 Code contributed by Neil Faught

#  RPresence package

library(readxl)
library(RPresence)
library(ggplot2)

# get the data read in
# Data for detections should be a data frame with rows corresponding to sites
# and columns to visits.
# The usual 1=detected; 0=not detected; NA=not visited is used.

input.data<-read.csv(file.path("..","pronghorn.csv"), header=TRUE, as.is=TRUE, strip.white=TRUE)
head(input.data)

input.history <- input.data[,c(2,3)] # the history extracted
head(input.history)

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.history)
ncol(input.history)
range(input.history, na.rm=TRUE) # check that all values are either 0 or 1
sum(is.na(input.history))    # are there any missing values?

# Get the sagebrush and aspect information. These are the covariates that will
# be used in the analysis. These are unit-level covariates.
unit.cov = input.data[,c(4,7)]
head(unit.cov)

# Create the *.pao file
prong.pao <- RPresence::createPao(input.history,
                                  unitcov=unit.cov,
                                  title='Pronghorn SSSS')
prong.pao

# define the list of models to fit
# Notice the commas between the column and the placement of the quotes
model.list.csv <- textConnection("
                                 p,               psi
                                 ~1,              ~1
                                 ~SURVEY,         ~1
                                 ~1,              ~aspect
                                 ~1,              ~sagebrush")

model.list <- read.csv(model.list.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)
model.list

# fit the models
model.fits <- plyr::alply(model.list, 1, function(x,detect.pao){
  cat("\n\n***** Starting ", unlist(x), "\n")
  fit <- RPresence::occMod(model=list(as.formula(paste("psi",x$psi)),
                                      as.formula(paste("p"  ,x$p  ))),
                           data=detect.pao,type="so")
  fit
},detect.pao=prong.pao)

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

# plot occupancy as a function of sagebrush density
psi.ma <- RPresence::modAvg(aic.table, param="psi")
head(psi.ma)
psi.ma$Site <- as.numeric(substring(row.names(psi.ma), 4+regexpr("unit",row.names(psi.ma), fixed=TRUE)))
plotdata <- data.frame(psi.ma, sagebrush = unit.cov$sagebrush)
head(plotdata)
ggplot(data=plotdata, aes(x=sagebrush, y=est))+
  ggtitle("Occupancy as a function of sagebrush density")+
  geom_point()+
  geom_ribbon(aes(ymin=lower_0.95, ymax=upper_0.95), alpha=0.2)+
  ylim(0,1)+
  ylab("Estimated occupancy")

# Plot occupancy as a function of aspect
plotdata <- data.frame(psi.ma, aspect = unit.cov$aspect)
head(plotdata)
ggplot(data=plotdata, aes(x=aspect, y=est))+
  ggtitle("Occupancy as a function of aspect")+
  geom_point()+
  geom_errorbar(aes(ymin=lower_0.95, ymax=upper_0.95), alpha=0.2)+
  ylim(0,1)+
  ylab("Estimated occupancy")
