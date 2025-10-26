# Single Species, Multi Season Occupancy analyais

#   House Finch example (ships with Presence)

# 2019-01-16 Code contributed by Neil Faught

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

input.data <- readxl::read_excel(file.path("..","HouseFinch.xlsx"), 
                                 sheet="Detections",
                                 na="-",
                                 col_names=FALSE)
head(input.data)

input.history <- input.data[,-1]
head(input.history)

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.history)
ncol(input.history)
range(input.history, na.rm=TRUE)
sum(is.na(input.history))

# Read in distance from Long Island (Site level covariate)
d <- readxl::read_excel(file.path("..","HouseFinch.xlsx"), 
                                 sheet="d",
                                 na="-",
                                 col_names=FALSE)
head(d)
site.cov <- d[,2]
names(site.cov) = 'd'

# Read in variable indicating whether house finches were detected on >10 stops in the previous year
f <- readxl::read_excel(file.path("..","HouseFinch.xlsx"), 
                                     sheet="f",
                                     na="-",
                                     col_names=FALSE)
head(f)
f <- f[,-1]
# Currently, this covariate is formatted as a survey level covariate when really we just need one
# entry per site, per season (as opposed to one entry per site, per visit, which is how the
# covariate is currently formatted).

# Only need to keep the 1st column out of every 50
#f <- f[,c(1,51,101,151,201,251)]
#head(f)
#names(f) = c("f1", "f2", "f3", "f4", "f5", "f6")

# Set all NAs to 0
#f[is.na(f)] = 0


# Create season covariate table (RPresence will technically treat this as a survey covariate).
# You need to stack the data by columns
season.cov <- data.frame(Site=1:nrow(f),
                         Visit=as.factor(rep(1:ncol(f), each=nrow(f))),
                         f =unlist(f),
                         stringsAsFactors=FALSE)

head(season.cov)

# Six seasons with 50 visits. Don't need same number of visits/season
Nvisits.per.season  <- rep(50,6) 

# Create the *.pao file
housefinch.pao <- RPresence::createPao(input.history,
                                    unitcov=site.cov,
                                    survcov=season.cov,
                                    nsurveyseason=Nvisits.per.season,
                                    title='housefinch SSMS')
housefinch.pao

# Define the models.
#    model.type do.1 is dynamic occupancy first parameterization
#               do.4 is dynamic occupancy 4th   parameterization (random occupancy)
# Random occupancy are fit using type="do.4" in the call.
#     Parameters are psi,  p with gamma=1-epsilon enforced internally

model.list.csv <- textConnection("
                                 p,               psi,   gamma,        epsilon,       model.type
                                 ~SEASON*d + f,   ~d,    ~SEASON*d,    ~d,            do.1
                                 ~SEASON*d + f,   ~d,    ~SEASON*d,    ~SEASON + d,   do.1
                                 ~SEASON*d + f,   ~d,    ~SEASON*d,    ~SEASON,       do.1
                                 ~SEASON*d + f,   ~d,    ~SEASON*d,    ~SEASON*d,     do.1
                                 ~SEASON*d + f,   ~d,    ~SEASON*d,    ~1,            do.1
                                 ~SEASON*d + f,   ~d,    ~d,           ~SEASON*d,     do.1
                                 ~SEASON*d + f,   ~d,    ~SEASON,      ~SEASON*d,     do.1
                                 ~SEASON*d + f,   ~d,    ~1,           ~SEASON*d,     do.1")


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
  
},detect.pao=housefinch.pao)


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


 