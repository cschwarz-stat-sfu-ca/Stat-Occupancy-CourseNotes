# Multi state single season

# Nichols et al (2007) Ecology 88:1395-1400 looked at breeding and
# non-breeding California Spotted Owls in s = 54 sites over k = 5
# surveys.

# Hand-fed mice conrmed occupancy; breeding status only
# confrmed when owl took mouse to nest.
#  . = site not visited.
#  0 = owl not detected.
#  1 = owl detected, but no detection of young.
#  2 = owl detected, along with evidence of young.

# Territories that were expected to be occupied were selectively monitored
# and so psi is expected to be close to 1 and not really of interest.

# See MacKenzie et al, Section 5.7.1 for more details

# 2018-10-01 Code contributed by Carl James Schwarz (cschwarz.stat.sfu.cs@gmail.com)

# General model averaging framework
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------

library(ggplot2)
library(plyr)
library(readxl)
library(reshape2)
library(RPresence)

# get the data
occ.data <- readxl::read_excel(file.path("..","CalSpottedOwl.xls"),
                               sheet="RawData",
                               skip=10, na='.')
head(occ.data)


input.history <- occ.data[,-1]
input.history[1:5,]

Nsites <- nrow(input.history)
Nvisits<- ncol(input.history)

# For biological reasons, we believe that errors in identifying states
# varies between the first 2 visits and the last 3 visit (early bs late
# breeding season). We create a site-time covariates.
# The format is the "long" format ordered by visit major order,i.e.
# stack the columns of the matrix of period covariates
#
#   The covariate matrix would be arranged originally as
#         1 2 3 4 5
#  site1 [E E L L L]
#  site2 [E E L L L]
#   :    [: : : : :]
#  siteN [E E L L L]
#
# and this is stacked by columns

period <- data.frame(Site=rep(1:Nsites, Nvisits),
                     Visit=rep(1:Nvisits, each=Nsites),
                     Per  =rep(c("E","E","L","L","L"), each=Nsites))
period[1:10,]

# create the pao file

owl.pao <- createPao(input.history,
                    survcov=period,
                    title="owl multi state single season")
summary(owl.pao)

#-----
# Use the multi-state model, but only for a single season 

# define the list of models to fit
# Notice the commas between the column and the placement of the quotes


#  psi	occupancy probability,
# r	 probability of being in the second state, conditional on the unit being occupied. 
# p	 detection probability. .
# delta	probability of detecting the second state in a survey, conditional on the species being detected in the survey. 
model.list.csv <- textConnection("
psi,     r,        p,     delta
~1,     ~1,       ~1,      ~1
~1,     ~1,       ~SURVEY, ~1
~1,     ~1,       ~STATE,  ~1
~1,     ~1,       ~1,      ~Per
~1,     ~1,       ~SURVEY,~Per
~1,     ~1,       ~STATE,  ~Per")

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
},detect.pao=owl.pao)




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



RPresence::modAvg(aic.table, param="psi")[1,]
RPresence::modAvg(aic.table, param="r")  [1,]

RPresence::modAvg(aic.table, param="p1")[seq(1,by=Nsites, length.out=Nvisits),]
RPresence::modAvg(aic.table, param="p2")[seq(1,by=Nsites, length.out=Nvisits),]
RPresence::modAvg(aic.table, param="delta")[seq(1,by=Nsites, length.out=Nvisits),]


