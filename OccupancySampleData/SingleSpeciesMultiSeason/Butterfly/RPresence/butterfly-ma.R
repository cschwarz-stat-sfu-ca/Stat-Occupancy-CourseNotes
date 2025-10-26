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

model.list.csv <- textConnection("
p,               psi,          gamma,      epsilon,  model.type
~1,              ~1,              ~1,           ~1,       do.1
~SEASON,         ~1,              ~Treatment,   ~Treatment, do.1
~SEASON,         ~1,              ~SEASON,      ~SEASON,  do.1")



model.list <- read.csv(model.list.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)
model.list


# fit the model
model.fits <- plyr::alply(model.list, 1, function(x,detect.pao){
  cat("\n\n***** Starting ", unlist(x), "\n")
  if(x$model.type == 'do.1'){
    fit <- RPresence::occMod(model=list(as.formula(paste("psi",x$psi)),
                                        as.formula(paste("p"  ,x$p  )),
                                        as.formula(paste("gamma",x$gamma)),
                                        as.formula(paste("epsilon",x$epsilon))),
                             data=detect.pao,type="do.1")
  }
  if(x$model.type == 'do.4'){
    fit <- RPresence::occMod(model=list(as.formula(paste("psi",x$psi)),
                                        as.formula(paste("p"  ,x$p  ))),
                             data=detect.pao,type="do.4")
  }
  fit <- RPresence.add.derived(fit)
  fit
  
},detect.pao=bfly.pao)


# Look at output from a specified model
model.number <- 3


names(model.fits[[model.number]])
names(model.fits[[model.number]]$real)
model.fits[[model.number]]$beta
names(model.fits[[model.number]]$derived)
model.fits[[model.number]]$derived$psi[1:10,]
model.fits[[model.number]]$real$gamma[1:5,]
model.fits[[model.number]]$real$epsilon[1:5,]

# Estimate of initial occupance
model.fits[[model.number]]$real$psi[grepl('unit1_', row.names(model.fits[[model.number]]$real$psi)),]

# Derived parameters - estimated occupancy for each unit in years 2....
names(model.fits[[model.number]]$derived)
model.fits[[model.number]]$derived$psi[ grepl('unit1_', row.names(model.fits[[model.number]]$derived$psi)),]

# Derived parameters - all of the psi stacked together
model.fits[[model.number]]$derived$all_psi[ grepl('unit1_', row.names(model.fits[[model.number]]$derived$all_psi)),]

# Estimate of  local extinction probability for each unit
model.fits[[model.number]]$real$epsilon[ seq(1, by=nrow(input.history), length.out=length(Nvisits.per.season)-1),]

# Estimate of  local colonization probability for each unit
model.fits[[model.number]]$real$gamma[ seq(1, by=nrow(input.history), length.out=length(Nvisits.per.season)-1),]

# Estimate of probability of detection at each time point for each unit
model.fits[[model.number]]$real$p[ grepl('unit1_', row.names(model.fits[[model.number]]$real$p), fixed=TRUE),]


# Get the change in occupancy
# Not yet possible to estimate the se of these values. May have to use bootstrapping.
model.fits[[model.number]]$derived$lambda [grepl('unit1_', row.names(model.fits[[model.number]]$derived$lambda),  fixed=TRUE),]
model.fits[[model.number]]$derived$lambdap[grepl('unit1_', row.names(model.fits[[model.number]]$derived$lambdap), fixed=TRUE),]



# collect models and make AIC table
aic.table <- RPresence::createAicTable(model.fits)
aic.table$table


# model averaging in the usual way

# initial occupancy
RPresence::modAvg(aic.table, param="psi")[1:5,]

# model averaging of derived parameters such as the occupancy at each time step
ma_all_psi <- RPresence.modAvg.derived(aic.table, param="all_psi")
psi.est <- ma_all_psi[grepl('unit1_', row.names(ma_all_psi),  fixed=TRUE),]
psi.est$Year <- as.numeric(substring(row.names(psi.est),1+regexpr("_",row.names(psi.est))))
psi.est$parameter <- 'psi'
psi.est


# likely more interested in colonization and extinction probabilities
epsilon.ma <- RPresence::modAvg(aic.table, param="epsilon")
epsilon.ma <- epsilon.ma[grepl('unit1$', row.names(epsilon.ma)),]
epsilon.ma$Year <- as.numeric(substr(row.names(epsilon.ma),7+regexpr("epsilon",row.names(epsilon.ma)),-1+regexpr("_",row.names(epsilon.ma))))
epsilon.ma$parameter <- 'epsilon'
epsilon.ma

gamma.ma <-RPresence::modAvg(aic.table, param="gamma")
gamma.ma <- gamma.ma[grepl('unit1$', row.names(gamma.ma)),]
gamma.ma$Year <- as.numeric(substr(row.names(gamma.ma),5+regexpr("gamma",row.names(gamma.ma)),-1+regexpr("_",row.names(gamma.ma))))
gamma.ma$parameter <- 'gamma'
gamma.ma

all.est <- rbind(psi.est, epsilon.ma, gamma.ma)
all.est$upper_0.95 = ifelse(all.est$upper_0.95 >1,1,all.est$upper_0.95)

ggplot(data=all.est, aes(x=Year,y=est, color=parameter))+
  ggtitle("Estimated occupancy, extinction, colonization, over time")+
  geom_point(position=position_dodge(w=0.2))+
  geom_line(position=position_dodge(w=0.2))+
  ylim(0,1)+
  geom_errorbar(aes(ymin=lower_0.95, ymax=upper_0.95), width=.1,position=position_dodge(w=0.2))+
  scale_x_continuous(breaks=1:3)+
  xlab("Season")

