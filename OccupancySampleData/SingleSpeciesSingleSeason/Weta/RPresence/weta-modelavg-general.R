# Single Species Single Season Occupancy 

# Weta Example with model averaging with a list of models to fit


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

input.history <- readxl::read_excel(file.path("..","weta.xls"),
                                    sheet="detection_histories",
                                    na="-",
                                    col_names=FALSE)  # notice no column names in row 1 of data file. 

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.history)
ncol(input.history)
range(input.history, na.rm=TRUE)
sum(is.na(input.history))

head(input.history)

# Get the site level covariates
site_covar <- readxl::read_excel(file.path("..","weta.xls"),
                                 sheet="site_covar",
                                 na="-",
                                 col_names=TRUE)  # notice col_names in row 1 of table. 


# Create an alternate site level covariate that is a categorical variable rather 
# than indicator variables
site_covar$BrowCat <- site_covar$BrowseCat
xtabs(~BrowCat, data=site_covar,exclude=NULL, na.action=na.pass)

head(site_covar)


# Get the individual covariates. 
obs1 <- readxl::read_excel(file.path("..","weta.xls"),
                           sheet="Obs1",
                           na="-",
                           col_names=FALSE) 
obs2 <- readxl::read_excel(file.path("..","weta.xls"),
                           sheet="Obs2",
                           na="-",
                           col_names=FALSE) 
obs3 <- readxl::read_excel(file.path("..","weta.xls"),
                           sheet="Obs3",
                           na="-",
                           col_names=FALSE) 

Obs <- obs1*1 + obs2*2 + obs3*3
head(Obs)

# Observational covariate needs to be "stacked" so that sites1...siteS for survey occastion 1
# are then followed by covariate at survey occastion 2 for sites1...siteS, etc

survey.cov <- data.frame(site=rep(1:nrow(input.history) , ncol(input.history)),
                         visit=as.character(rep(1:ncol(input.history), each=nrow(input.history))),  # notice we make a character 
                         obs1 =as.vector(unlist(obs1)),
                         obs2 =as.vector(unlist(obs2)),
                         obs3 =as.vector(unlist(obs3)),
                         Obs  =paste("O",as.vector(unlist(Obs)),sep=""),    # notice we make a character string
                         stringsAsFactors=FALSE)
survey.cov$Obs[ grepl("NA",survey.cov$Obs)] <- NA
head(survey.cov[,c("visit","site","Obs")],n=100)

str(survey.cov)

# check that missing values in history and observer covariates align
select <- is.na(as.vector(unlist(input.history)))
survey.cov[select,]
sum(is.na(survey.cov[!select,]))


# Create the *.pao file
weta.pao <- RPresence::createPao(input.history,
                                 unitcov=site_covar,
                                 survcov=survey.cov,
                                 title='weta SSSS')
weta.pao


# define the list of models to fit
# Notice the commas between the column and the placement of the quotes
model.list.csv <- textConnection("
p,               psi
~1,              ~1
~1,              ~BrowCat
~BrowCat,        ~BrowCat
~visit,          ~BrowCat
~Obs,            ~BrowCat
~Obs+visit,      ~BrowCat
~Obs+visit+BrowCat,      ~BrowCat
~Obs+visit+BrowCat,      ~1
~Obs+visit,      ~1")


model.list <- read.csv(model.list.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)
model.list


# fit the models
model.fits <- plyr::alply(model.list, 1, function(x,detect.pao){
  cat("\n\n***** Starting ", unlist(x), "\n")
  fit <- RPresence::occMod(model=list(as.formula(paste("psi",x$psi)),
                                      as.formula(paste("p"  ,x$p  ))),
                           data=detect.pao,type="so")
  fit
},detect.pao=weta.pao)




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




