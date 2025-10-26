# multi season multi-state model

# Robust design, Multi State, Reproduction parameterization
# Using RPresence

# Visits to hawk territories to establish occupancy and breeding
# 6 years x up to 9 visits/year

# 2018-11-27 Code contributed by Neil Faught


#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------

library(car)
library(ggplot2)
library(readxl)
library(reshape2)
library(RPresence)

input.data <- read.csv(file.path("..","hawk.csv"),
                          as.is=TRUE, strip.white=TRUE, header=TRUE)
head(input.data)

input.history = input.data[,-c(1,2)]
head(input.history)

# Create site level covariate data frame
unit.cov = input.data[,c(1,2)]
names(unit.cov) = c("Territory","tsh")
head(unit.cov)

# make a plot of the imputed occupancy
plotdata <- reshape2::melt(input.data,
                           id.vars=c("TerritoryN","tsh"),
                           value.name="ObsOcc",
                           variable.name="Visit")
plotdata$Visit <- as.numeric(substring(plotdata$Visit,2))
head(plotdata)

ggplot(data=plotdata, aes(x=Visit, y=TerritoryN, color=as.factor(ObsOcc), shape=as.factor(ObsOcc)))+
  ggtitle("Observed occupancy state", subtitle="Not corrected for false negatives")+
  geom_point()+
  scale_shape_manual(name="Observed\nOccupancy", values=c(1, 16, 19))+
  scale_color_manual(name="Observed\nOccupancy", values=c("black","green","red"))+
  theme_bw()

# Number of visits in each season
Nvisits.per.season  <- rep(6,9) 

# create the pao file
hawk.pao <- createPao(input.history,
                      nsurveyseason=Nvisits.per.season,
                      unitcov=unit.cov,
                      title="hawk multi state multi season")
summary(hawk.pao)



#  psi	occupancy probability,
# r	 probability of being in the second state, conditional on the unit being occupied. 
# p	 detection probability. .
# delta	probability of detecting the second state in a survey, conditional on the species being detected in the survey. 
model.list.csv <- textConnection("
psi,         r,      p,          delta
~DYN + PREV_STATE,       ~DYN+ PREV_STATE,     ~1,         ~1
~tsh*DYN,   ~tsh,   ~1,         ~1
~tsh*DYN,   ~1,     ~1,         ~1")

model.list <- read.csv(model.list.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)
model.list


# fit the models
model.fits <- plyr::alply(model.list, 1, function(x,detect.pao){
  cat("\n\n***** Starting ", unlist(x), "\n")
  fit <- RPresence::occMod(model=list(as.formula(paste("psi",  x$psi)),
                                      as.formula(paste("r"    ,x$r  )),
                                      as.formula(paste("p"    ,x$p  )),
                                      as.formula(paste("delta",x$delta))),
                           data=detect.pao,type="do.ms.2")
  fit
},detect.pao=hawk.pao)

# Look the output from a specific model
check.model <- 1

names(model.fits[[check.model]])
model.fits[[check.model]]$beta

names(model.fits[[check.model]]$real)
model.fits[[check.model]]$real$psi[1:5,]
model.fits[[check.model]]$real$Cpsi0[1:5,]
model.fits[[check.model]]$real$Cpsi1[1:5,]
model.fits[[check.model]]$real$Cpsi2[1:5,]
model.fits[[check.model]]$real$r[1:5,]
model.fits[[check.model]]$real$CR0[1:5,]
model.fits[[check.model]]$real$CR1[1:5,]
model.fits[[check.model]]$real$CR2[1:5,]
model.fits[[check.model]]$real$delta[1:5,]

model.fits[[check.model]]$real$p1[1:5,]
model.fits[[check.model]]$real$p2[1:5,]

names(model.fits[[check.model]]$derived)

# Model averaging
aic.table <- RPresence::createAicTable(model.fits)
aic.table$table

names(aic.table)
