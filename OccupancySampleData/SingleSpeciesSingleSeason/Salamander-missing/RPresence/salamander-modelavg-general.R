# Single Species Single Season Occupancy 

# Salamander Example with model averaging with a list of models to fit

# 2018-08-15 Code contributed by Carl James Schwarz (cschwarz.stat.sfu.cs@gmail.com)

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
library(plyr)
library(readxl)
library(RPresence)
library(ggplot2)

# get the data read in
# Data for detections should be a data frame with rows corresponding to sites
# and columns to visits.
# The usual 1=detected; 0=not detected; NA=not visited is used.

input.data <- readxl::read_excel("../salamander.xls",
                                 sheet="MissingData",
                                 na="-",
                                 col_names=FALSE)  # notice no column names in row 1 of data file. 

head(input.data)


# Extract the history records
input.history <- input.data[, 1:5] # the history extracted
head(input.history)



# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.history)
ncol(input.history)
range(input.history, na.rm=TRUE) # check that all values are either 0 or 1
sum(is.na(input.history))    # are there any missing values?


# Create the survey covariates
# This covariate needs to be "stacked" so that sites1...site39 for survey occastion 1
# are then followed by covariate at survey occastion 2 for sites1...site39, etc

survey.cov <- data.frame(site=rep(1:nrow(input.history), ncol(input.history)),
                         visit=rep(1:ncol(input.history), each=nrow(input.history)),
                         d=rep( c("d1","d1","d2","d2","d2"),  each=nrow(input.history)))
head(survey.cov)

survey.cov[c(1:4, 37:41),]

# Create the *.pao file
salamander.pao <- RPresence::createPao(input.history,
                                       title='Salamander SSSS',
                                       survcov=survey.cov)
salamander.pao
summary(salamander.pao)



# define the list of models to fit
# Notice the commas between the column and the placement of the quotes
model.list.csv <- textConnection("
p,               psi
~1,              ~1
~SURVEY,        ~1
~d,             ~1")

model.list <- read.csv(model.list.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)
model.list


# fit the models
model.fits <- plyr::alply(model.list, 1, function(x,detect.pao){
  cat("\n\n***** Starting ", unlist(x), "\n")
  fit <- RPresence::occMod(model=list(as.formula(paste("psi",x$psi)),
                                      as.formula(paste("p"  ,x$p  ))),
                           data=detect.pao,type="so")
  fit
},detect.pao=salamander.pao)




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



RPresence::modAvg(aic.table, param="psi")[1:5,]


ma.p <- RPresence::modAvg(aic.table, param="p")
ma.p[grepl("unit1$", rownames(ma.p)),]




