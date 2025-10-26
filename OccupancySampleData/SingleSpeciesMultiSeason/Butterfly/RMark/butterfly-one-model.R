#Analysis of butterfly data set
#3 seasons, 2 sampling occasions per season
#Interested in occupancy of CLLESSU
#Using RMark

library(car)
library(readxl)
library(RMark)
library(ggplot2)
library(reshape)

# Get the RMark additional functions 
source(file.path("..","..","..","AdditionalFunctions","RMark.additional.functions.R"))

#Read in data
bfly<-read.csv("../Butterfly.csv")
#Data is in long format, need to convert to wide format for RMark
#Want one row for each transect
head(bfly)
bfly<-subset(bfly,select=c(Field,Border,Transect,Treatment,Date,Visit,CLLESSU))
head(bfly)
bfly$Transect=as.factor(bfly$Transect)

#Convert to wide format
junk<-melt(bfly,id.var=c("Transect","Visit"),measure.var="CLLESSU")
j2=cast(junk,Transect ~ Visit)
head(j2)
#Note how there are counts in these columns, need to convert to simply 0/1 for 
#present/absent
j2$`1`=ifelse(j2$`1`>0,1,0)
j2$`2`=ifelse(j2$`2`>0,1,0)
j2$`3`=ifelse(j2$`3`>0,1,0)
j2$`4`=ifelse(j2$`4`>0,1,0)
j2$`5`=ifelse(j2$`5`>0,1,0)
j2$`6`=ifelse(j2$`6`>0,1,0)
head(j2)

input.history = j2[,2:7]
# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.history)
ncol(input.history)
range(input.history, na.rm=TRUE)
sum(is.na(input.history))

#Format the capture history to be used by RMark
input.history <- data.frame(freq=1,
                            ch=apply(input.history,1,paste, collapse=""), stringsAsFactors=FALSE)
head(input.history)
# Change any NA to . in the chapter history
input.history$ch <- gsub("NA",".", input.history$ch, fixed=TRUE)
head(input.history)

#Is the treatment a site or survey level covariate (i.e. does each transect only
#have a single treatment applied to it over the course of the 6 sampling occasions?)
jcov = melt(bfly,id.var=c("Transect","Visit"),measure.var="Treatment")
jcov2=cast(jcov,Transect ~ Visit)
head(jcov2)
jcov2$test = ifelse(jcov2$`1`==jcov2$`2` & jcov2$`1`==jcov2$`3` &
                   jcov2$`1`==jcov2$`4` & jcov2$`1`==jcov2$`5` &
                   jcov2$`1`==jcov2$`6`,0,1)
head(jcov2)
#Looks like this is a site level covariate
sum(jcov2$test)

#Add site covariates to input history
input.history = cbind(input.history,jcov2$`1`)
names(input.history)[3] = "Treatment"
head(input.history)



#Create the RMark data structure
#2 Seasons, with 3 visits per season
max.visit.per.year <- 2
n.season <- 3

bfly.data <- process.data(data=input.history, 
                           model="RDOccupEG",
                           time.intervals=c( rep( c(rep(0,max.visit.per.year-1),1),n.season-1),
                                             rep(0,max.visit.per.year-1)),
                          group = "Treatment")
summary(bfly.data)

# add visit level covariates
# In this case none
bfly.ddl <- make.design.data(bfly.data)
bfly.ddl

#Fit Null Model
mod.fit <-  RMark::mark(bfly.data,
                        ddl = bfly.ddl,
                        model="RDOccupEG",
                        model.parameters=list(
                          Psi   =list(formula=~1),
                          p     =list(formula=~1),
                          Epsilon = list(formula=~1),
                          Gamma = list(formula=~1)
                        )
)

# Estimates of real initial occupancy, detection, local exctinction,
# and local colonization probabilities
mod.fit$results$real
# Estimates of regression coefficients
mod.fit$results$beta
# Estimates of various derived parameters for estimated occpancy in each year,
# estimated ratio of in occupancy from year to year, and estimated log
# of the ratio of odds in occupancy from year to year
mod.fit$results$derived

trt = c(rep("burn",3),rep("burncontrol",3),rep("disk",3),rep("diskcontrol",3),rep("fieldcontrol",3))
trt2 = c(rep("burn",2),rep("burncontrol",2),rep("disk",2),rep("diskcontrol",2),rep("fieldcontrol",2))

mod.fit$results$derived$"psi Probability Occupied"$Treatment = trt
mod.fit$results$derived$"psi Probability Occupied"
mod.fit$results$derived$"lambda Rate of Change"$Treatment = trt2
mod.fit$results$derived$"lambda Rate of Change"
mod.fit$results$derived$"log odds lambda"$Treatment = trt2
mod.fit$results$derived$"log odds lambda"

# Plot derived occupancy probabilities
plotdata = mod.fit$results$derived$"psi Probability Occupied"
plotdata$Season = rep(c(1:3),5)
head(plotdata)


ggplot(data=plotdata, aes(x=Season,y=estimate))+
  ggtitle("Estimated occupancy over time")+
  geom_point(position=position_dodge(w=0.2))+
  geom_line(position=position_dodge(w=0.2))+
  ylim(0,1)+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.1,position=position_dodge(w=0.2))+
  scale_x_continuous(breaks=1:3)+
  facet_wrap(~Treatment)

cleanup(ask=FALSE)
