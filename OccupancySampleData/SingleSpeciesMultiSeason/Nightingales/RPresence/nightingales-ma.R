# Single Species, Multi Season Occupancy analyais

#   Nighingale detection..
#   55 sites x 10 years x 8 visits/site
#   No covariates

# 2018-11-15 Code contributed by Carl James Schwarz (cschwarz.stat.sfu.cs@gmail.com)

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

#  RPresence package

library(readxl)
library(RPresence)
library(ggplot2)

# Get the RPResence additional functions 
source(file.path("..","..","..","AdditionalFunctions","Rpresence.additional.functions.R"))

# get the data read in
# Data for detections should be a data frame with rows corresponding to sites
# and columns to visits.
# The usual 1=detected; 0=not detected; NA=not visited is used.

input.data <- read.csv(file.path("..","nightingales.csv"), header=TRUE, as.is=TRUE, strip.white=TRUE)
head(input.data)



input.history <- input.data[, -1]
head(input.history)

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.history)
ncol(input.history)
range(input.history, na.rm=TRUE)
sum(is.na(input.history))

# site covariates - none.
#site.covar <- NULL





# Number of visits in each season
Nvisits.per.season  <- rep(8,10) 

# Create the *.pao file
nightingale.pao <- RPresence::createPao(input.history,
                                nsurveyseason=Nvisits.per.season,
                                #unitcov=site.covar,  # no site covariates
                                title='Nighingale SSMS')
nightingale.pao
summary(nighingale.pao)


# Define the models.
#    model.type do.1 is dynamic occupancy first parameterization
#               do.4 is dynamic occupancy 4th   parameterization (random occupancy)
# Random occupancy are fit using type="do.4" in the call.
#     Parameters are psi,  p with gamma=1-epsilon enforced internally

model.list.csv <- textConnection("
p,               psi,          gamma,      epsilon,  model.type
~1,              ~1,              ~1,           ~1,       do.1
~SEASON,         ~1,              ~1,           ~1,       do.1
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

},detect.pao=nightingale.pao)


# Look at output from a specified model
model.number <- 2


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

ggplot(data=all.est, aes(x=Year,y=est, color=parameter))+
   ggtitle("Estimated occupancy, extinction, colonization, over time")+
   geom_point(position=position_dodge(w=0.2))+
   geom_line(position=position_dodge(w=0.2))+
   ylim(0,1)+
   geom_errorbar(aes(ymin=lower_0.95, ymax=upper_0.95), width=.1,position=position_dodge(w=0.2))+
   scale_x_continuous(breaks=1:10)


