#Analysis of butterfly data set
#3 seasons, 2 sampling occasions per season
#Interested in occupancy of CLLESSU
#Using RPresence
#2018-11-26 submitted by Neil Faught

library(car)
library(readxl)
library(RPresence)
library(ggplot2)
library(reshape)

# Get the RMark additional functions 
source(file.path("..","..","..","AdditionalFunctions","RPresence.additional.functions.R"))

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
site.covar = data.frame(Treatment = jcov2$`1`)
head(site.covar)

# Number of visits in each season
Nvisits.per.season  <- rep(2,3) 

# Create the *.pao file
bfly.pao <- RPresence::createPao(input.history,
                                       nsurveyseason=Nvisits.per.season,
                                       unitcov=site.covar,
                                       title='Nighingale SSMS')
bfly.pao
summary(bfly.pao)

#Fit Null Model
mod.fit <- RPresence::occMod(
  model=list(psi~1, gamma~1, epsilon~1, p~1), 
  type="do.1", data=bfly.pao)
summary(mod.fit)

mod.fit <- RPresence.add.derived(mod.fit)


# Estimate of initial occupance
mod.fit$real$psi[1,]

# Derived parameters - estimated occupancy for each unit in years 2....
names(mod.fit$derived)
mod.fit$derived$psi[ grepl('unit1_', row.names(mod.fit$derived$psi)),]

# Additional derived parameters - all of the psi stacked together
mod.fit$derived$all_psi[ grepl('unit1_', row.names(mod.fit$derived$all_psi)),]

# Estimate of  local extinction probability for each unit
mod.fit$real$epsilon[ seq(1, by=nrow(input.history), length.out=length(Nvisits.per.season)-1),]

# Estimate of  local colonization probability for each unit
mod.fit$real$gamma[ seq(1, by=nrow(input.history), length.out=length(Nvisits.per.season)-1),]

# Estimate of probability of detection at each time point for each unit
mod.fit$real$p[ grepl('_unit1$', row.names(mod.fit$real$p)),]

# Derived parameters - estimated occupancy for each unit in years 2....
names(mod.fit$derived)
mod.fit$derived$lambda[ grepl('unit1_', row.names(mod.fit$derived$lambda)),]


# Get the change in occupancy
# Not yet possible to estimate the se of these values. May have to use bootstrapping.
mod.fit$derived$lambda [grepl('unit1_', row.names(mod.fit$derived$lambda),  fixed=TRUE),]
mod.fit$derived$lambdap[grepl('unit1_', row.names(mod.fit$derived$lambdap), fixed=TRUE),]

# Plot each season's predicted occupancy
plotdata = mod.fit$derived$all_psi[ grepl('unit1_', row.names(mod.fit$derived$all_psi)),]
plotdata$Season = c(1:3)
head(plotdata)

ggplot(data=plotdata, aes(x=Season,y=est))+
  ggtitle("Estimated occupancy over time")+
  geom_point(position=position_dodge(w=0.2))+
  geom_line(position=position_dodge(w=0.2))+
  ylim(0,1)+
  geom_errorbar(aes(ymin=lower_0.95, ymax=upper_0.95), width=.1,position=position_dodge(w=0.2))+
  scale_x_continuous(breaks=1:3)
