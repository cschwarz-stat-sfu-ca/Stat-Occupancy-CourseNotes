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

#
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
#
#  psi	occupancy probability,
# r	 probability of being in the second state, conditional on the unit being occupied. 
# p	 detection probability. .
# delta	probability of detecting the second state in a survey, conditional on the species being detected in the survey. 


mod1 <- occMod(model=list(
    psi~1, # occupancy regardless of state
    r~1,   # occupancy in state 2 | occupied
    p~1,   # detection in state i
    delta~1),# identified as state 2 if detected (and in state 2)
    data=owl.pao,
    type="do.ms.2")
summary(mod1)

names(mod1$real)

mod1$real$psi[1,] # prob site is occupied
mod1$real$r[1,]  # prob in state 2 | occupied


mod1$real$p1[seq(1,by=Nsites, length.out=Nvisits),]  # prob detect if in state 1
mod1$real$p2[seq(1,by=Nsites, length.out=Nvisits),]  # prob detect if in state 2
mod1$real$delta[seq(1,by=Nsites, length.out=Nvisits),] # prob identify in in state 2 if detect in state 2

#-------------------------------------------------------------
# Need a sutiable model for the detection probabilities


mod2 <- occMod(model=list(
    psi~1, # occupancy regardless of state
    r~1,   # occupancy in state 2 | occupied
    p~SURVEY,   # detection in each state vary by survey
    delta~1),# identified as state 2 if detected (and in state 2)
    data=owl.pao,
    type="do.ms.2")
summary(mod2)

names(mod2$real)

mod2$real$psi[1,] # prob site is occupied
mod2$real$r[1,]   # prob in state 2 | occupied


mod2$real$p1[seq(1,by=Nsites, length.out=Nvisits),]  # prob detect if in state 1
mod2$real$p2[seq(1,by=Nsites, length.out=Nvisits),]  # prob detect if in state 2
mod2$real$delta[seq(1,by=Nsites, length.out=Nvisits),] # prob identify in in state 2 if detect in state 2
 
#-------------------------------------------------------------
# Need a sutiable model for the detection probabilities


mod3 <- occMod(model=list(
    psi~1, # occupancy regardless of state
    r~1,   # occupancy in state 2 | occupied
    p~STATE,   # detection in each state vary by state
    delta~1),# identified as state 2 if detected (and in state 2)
    data=owl.pao,
    type="do.ms.2")
summary(mod3)

names(mod3$real)

mod3$real$psi[1,] # prob site is occupied
mod3$real$r[1,]   # prob in state 2 | occupied


mod3$real$p1[seq(1,by=Nsites, length.out=Nvisits),]  # prob detect if in state 1
mod3$real$p2[seq(1,by=Nsites, length.out=Nvisits),]  # prob detect if in state 2
mod3$real$delta[seq(1,by=Nsites, length.out=Nvisits),] # prob identify in in state 2 if detect in state 2
 

#-------------------------------------------------------------
# Need a sutiable model for the detection probabilities

# On the surface, this model does not appear to converge, but the problem is
# tht delta is estimated to be 0 in the first period which causes numerical
# problems (logit(0) is -infinity.)

mod4 <- occMod(model=list(
    psi~1, # occupancy regardless of state
    r~1,   # occupancy in state 2 | occupied
    p~1,   # detection in each state vary by state
    delta~Per),# identified as state 2 if detected (and in state 2)
    data=owl.pao,
    type="do.ms.2")
summary(mod4)

names(mod4$real)
mod4$beta$delta

mod4$real$psi[1,] # prob site is occupied
mod4$real$r[1,]   # prob in state 2 | occupied


mod4$real$p1[seq(1,by=Nsites, length.out=Nvisits),]  # prob detect if in state 1
mod4$real$p2[seq(1,by=Nsites, length.out=Nvisits),]  # prob detect if in state 2
mod4$real$delta[seq(1,by=Nsites, length.out=Nvisits),] # prob identify in in state 2 if detect in state 2
 
#-------------------------------------------------------------
# Need a sutiable model for the detection probabilities

mod5 <- occMod(model=list(
    psi~1, # occupancy regardless of state
    r~1,   # occupancy in state 2 | occupied
    p~SURVEY,   # detection in each state vary by state
    delta~Per),# identified as state 2 if detected (and in state 2)
    data=owl.pao,
    type="do.ms.2")
summary(mod5)

names(mod5$real)

mod5$real$psi[1,] # prob site is occupied
mod5$real$r[1,]   # prob in state 2 | occupied


mod5$real$p1[seq(1,by=Nsites, length.out=Nvisits),]  # prob detect if in state 1
mod5$real$p2[seq(1,by=Nsites, length.out=Nvisits),]  # prob detect if in state 2
mod5$real$delta[seq(1,by=Nsites, length.out=Nvisits),] # prob identify in in state 2 if detect in state 2
 

#-------------------------------------------------------------
# Need a sutiable model for the detection probabilities

mod6 <- occMod(model=list(
    psi~1, # occupancy regardless of state
    r~1,   # occupancy in state 2 | occupied
    p~STATE,   # detection in each state vary by state
    delta~Per),# identified as state 2 if detected (and in state 2)
    data=owl.pao,
    type="do.ms.2")
summary(mod6)

names(mod6$real)

mod6$real$psi[1,] # prob site is occupied
mod6$real$r[1,]   # prob in state 2 | occupied


mod6$real$p1[seq(1,by=Nsites, length.out=Nvisits),]  # prob detect if in state 1
mod6$real$p2[seq(1,by=Nsites, length.out=Nvisits),]  # prob detect if in state 2
mod6$real$delta[seq(1,by=Nsites, length.out=Nvisits),] # prob identify in in state 2 if detect in state 2
 


#-------------
# Model averaging
models<-list(mod1, mod2, mod3, 
             mod4, mod5, mod6 
            )
results<-RPresence::createAicTable(models)
summary(results)

RPresence::modAvg(results, param="psi")[1,]
RPresence::modAvg(results, param="r")  [1,]

RPresence::modAvg(results, param="p1")[seq(1,by=Nsites, length.out=Nvisits),]
RPresence::modAvg(results, param="p2")[seq(1,by=Nsites, length.out=Nvisits),]
RPresence::modAvg(results, param="delta")[seq(1,by=Nsites, length.out=Nvisits),]

