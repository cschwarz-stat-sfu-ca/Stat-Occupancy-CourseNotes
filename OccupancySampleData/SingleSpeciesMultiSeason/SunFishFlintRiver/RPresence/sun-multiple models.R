# Sunfish
# Redbreast sunfish occupancy of streams segments over 8 seasons = 4 years * summer*spring
# 3 quadrats were sampled per segment via electrofishing total 24 sampling occasions
# 15 covariates: LN(link magnitude) and standardized Q (streamflow) for the sampling intervals;
# minimum discharge for intervals 1-7, maximum discharge intervals 1-7
  
# Fitting several models using the RPresence Package
# Single Species Multiple Seasons

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

library(readxl)
library(RPresence)
library(ggplot2)

# Get the RPResence additional functions 
source(file.path("..","..","..","AdditionalFunctions","Rpresence.additional.functions.R"))

# get the data read in
# Data for detections should be a data frame with rows corresponding to sites
# and columns to visits.
# The usual 1=detected; 0=not detected; NA=not visited is used.

input.data <- readxl::read_excel(file.path("..","sunfish.xls"), 
                                 skip=1, col_names=TRUE)#No NAs in this data set, make sure to use "na =" statement if there are
head(input.data)

input.history <- input.data[, 1:24]
head(input.history)
# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.history)
ncol(input.history)
range(input.history, na.rm=TRUE)
sum(is.na(input.history))

#Create data frame of site covariates
site.cov <- data.frame(Site=1:nrow(input.data), input.data[,25])
head(site.cov)

#Create data frame of survey covariates
minQ.data = input.data[,26:32]
minQ.data = cbind(minQ.data,rep(NA,nrow(minQ.data)))
names(minQ.data)[8] = "minQ8"
head(minQ.data)
maxQ.data = input.data[,33:39]
maxQ.data = cbind(maxQ.data,rep(NA,nrow(maxQ.data)))
names(maxQ.data)[8] = "maxQ8"
head(maxQ.data)

# Create survey covariate table. You need to stack the data by columns
survey.cov <- data.frame(Site=1:nrow(minQ.data),
                         Visit=as.factor(rep(1:ncol(minQ.data), each=nrow(minQ.data))),
                         minQ =unlist(minQ.data),
                         maxQ =unlist(maxQ.data),
                         stringsAsFactors=FALSE)

head(survey.cov)

# Eight seasons with 3 visits. Don't need same number of visits/season
Nvisits.per.season  <- rep(3,8) # 8 seasons with 3 visits. Don't need same number of visits/season

# Create the *.pao file
sunfish.pao <- RPresence::createPao(input.history,
                                  unitcov=site.cov,
                                  survcov=survey.cov,
                                  nsurveyseason=Nvisits.per.season,
                                  title='sunfish SSMS')
sunfish.pao

# Define the models.
#    model.type do.1 is dynamic occupancy first parameterization
#               do.4 is dynamic occupancy 4th   parameterization (random occupancy)
# Random occupancy are fit using type="do.4" in the call.
#     Parameters are psi,  p with gamma=1-epsilon enforced internally

model.list.csv <- textConnection("
                                 p,               psi,         gamma,      epsilon,  model.type
                                 ~SEASON,         ~LN,         ~maxQ,      ~minQ,       do.1
                                 ~SEASON,         ~LN,         ~minQ,      ~minQ,       do.1
                                 ~SEASON,         ~LN,         ~minQ,      ~maxQ,       do.1
                                 ~SEASON,         ~1,         ~SEASON,    ~SEASON,     do.1")


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
  
},detect.pao=sunfish.pao)

# Look at output from a specified model
model.number <- 1


names(model.fits[[model.number]])
names(model.fits[[model.number]]$real)
model.fits[[model.number]]$beta
names(model.fits[[model.number]]$derived)
model.fits[[model.number]]$derived$psi[1:10,]
model.fits[[model.number]]$real$gamma[1:5,]
model.fits[[model.number]]$real$epsilon[1:5,]
# Estimate of initial occupance
model.fits[[model.number]]$real$psi[grepl('unit1_', row.names(model.fits[[model.number]]$real$psi)),]
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
