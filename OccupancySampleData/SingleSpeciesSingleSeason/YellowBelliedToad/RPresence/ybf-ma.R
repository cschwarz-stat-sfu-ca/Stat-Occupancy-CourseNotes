# Single Species Single Season Occupancy 

# Yellow-bellied toad
# Single Season Single Season occupancy 

#  RPresence package

# 2018-12-02 Code contributed by Neil Faught

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------



library(readxl)
library(RPresence)
library(ggplot2)

# get the data read in
# Data for detections should be a data frame with rows corresponding to sites
# and columns to visits.
# The usual 1=detected; 0=not detected; NA=not visited is used.

input.history <- readxl::read_excel(file.path("..","YellowBelliedToad.xlsx"),
                                    sheet="detections")

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.history)
ncol(input.history)
range(input.history, na.rm=TRUE)
sum(is.na(input.history))

head(input.history)

# Get the site level covariates - none.



# Get the site x visit (sampling) covariates. 
jdate <- readxl::read_excel(file.path("..","YellowBelliedToad.xlsx"),
                            sheet="JulianDate")
range(jdate)

# Observational covariate needs to be "stacked" so that sites1...siteS for survey occastion 1
# are then followed by covariate at survey occastion 2 for sites1...siteS, etc

survey.cov <- data.frame(site =rep(1:nrow(input.history) , ncol(input.history)),
                         visit =as.character(rep(1:ncol(input.history), each=nrow(input.history))),  # notice we make a character 
                         jdate =as.vector(unlist(jdate)), 
                         stringsAsFactors=FALSE)
head(survey.cov[,c("visit","site","jdate")],n=50)
str(survey.cov)

# standarize jdate
survey.cov$jdateS <- (survey.cov$jdate - 150)/10

# check that missing values in history and observer covariates align
select <- is.na(as.vector(unlist(input.history)))
survey.cov[select,]
sum(is.na(survey.cov[!select,]))


# Create the *.pao file
ybf.pao <- RPresence::createPao(input.history,
                                #unitcov=site_covar, # no site covariates
                                survcov=survey.cov, 
                                title='ybf SSSS')
ybf.pao



# define the list of models to fit
# Notice the commas between the column and the placement of the quotes
model.list.csv <- textConnection("
                                 p,               psi
                                 ~1,              ~1
                                 ~SURVEY,         ~1
                                 ~jdateS,         ~1
                                 ~jdateS + I(jdateS^2),          ~1")

model.list <- read.csv(model.list.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)
model.list


# fit the models
model.fits <- plyr::alply(model.list, 1, function(x,detect.pao){
  cat("\n\n***** Starting ", unlist(x), "\n")
  fit <- RPresence::occMod(model=list(as.formula(paste("psi",x$psi)),
                                      as.formula(paste("p"  ,x$p  ))),
                           data=detect.pao,type="so")
  fit
},detect.pao=ybf.pao)

# Look the output from a specific model
check.model <- 4

names(model.fits[[check.model]])
model.fits[[check.model]]$beta

names(model.fits[[check.model]]$real)
model.fits[[check.model]]$real$psi[1:5,]
model.fits[[check.model]]$real$p[1:5,]

names(model.fits[[check.model]]$derived)
model.fits[[check.model]]$derived$psi_c[1:10,]
tail(model.fits[[check.model]]$derived$psi_c)

#------
# Model averaging
results<-RPresence::createAicTable(model.fits)
summary(results)

RPresence::modAvg(results, param="psi")[1:5,]

# look at detectability 
ma.p <- RPresence::modAvg(results, param="p")

# we need to extract the site and survey and merge with survey.cov to get jdate
head(ma.p)
ma.p$site <- as.numeric( substring(row.names(ma.p), 4+regexpr("unit", row.names(ma.p))))
ma.p$visit<- substr(row.names(ma.p),2,2)
head(ma.p)

ma.p <- merge(ma.p, survey.cov)
head(ma.p)

ggplot(data=ma.p, aes(x=jdate, y=est))+
  ggtitle("Estimated detection as a function of julian date")+
  geom_point()+
  geom_ribbon(aes(ymin=lower_0.95, ymax=upper_0.95), alpha=0.2)+
  ylim(0,1)

