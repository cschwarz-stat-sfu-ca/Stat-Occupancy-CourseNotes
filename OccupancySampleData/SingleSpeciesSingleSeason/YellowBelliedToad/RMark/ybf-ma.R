# Single Species Single Season Occupancy 

# Yellow-bellied toad
# Single Season Single Season occupancy 

#  RMark package

# 2018-12-02 Code contributed by Neil Faught

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------



library(readxl)
library(RMark)
library(ggplot2)

# get the data read in
# Data for detections should be a data frame with rows corresponding to sites
# and columns to visits.
# The usual 1=detected; 0=not detected; NA=not visited is used.

input.detect <- readxl::read_excel(file.path("..","YellowBelliedToad.xlsx"),
                                   sheet="detections")

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.detect)
ncol(input.detect)
range(input.detect, na.rm=TRUE)
sum(is.na(input.detect))

head(input.detect)


# Extract the history records and create a capture history
input.history <- data.frame(freq=1,
                            ch=apply(input.detect[,1:2],1,paste, collapse=""), stringsAsFactors=FALSE)
head(input.history)

# Change any NA to . in the chapter history
input.history$ch <- gsub("NA",".", input.history$ch, fixed=TRUE)
head(input.history)


# Get the site level covariates - none.



# Get the site x visit (sampling) covariates. 
jdate <- readxl::read_excel(file.path("..","YellowBelliedToad.xlsx"),
                            sheet="JulianDate")
range(jdate)

# standarize jdate
jdateS <- (jdate - 150)/10

names(jdate ) <- c('jdate1'   ,'jdate2')
names(jdateS) <- c('jdateS1'  ,'jdateS2')

# Sampling covariates must be added to the input.history as multiple columns
# We also need to create the quadratic terms in advance - this is a pain
jdateSQ <- jdateS
names(jdateSQ)<- c('jdateSQ1' , 'jdateSQ2')
jdateSQ <- jdateSQ^2


input.history <- cbind(input.history, jdate, jdateS, jdateSQ)
head(input.history)

ybf.data <- process.data(data=input.history,
                         model="Occupancy")
summary(ybf.data)

# Set up survey covariates in the ddl (in this case none)
ybf.ddl <- make.design.data(ybf.data)
ybf.ddl

# What are the parameter names for Single Season Single Species models
setup.parameters("Occupancy", check=TRUE)


# define the list of models to fit
# Notice the commas between the column and the placement of the quotes
model.list.csv <- textConnection("
                                 p,         Psi
                                 ~1,              ~1
                                 ~time,           ~1
                                 ~jdateS,         ~1
                                 ~jdateS+jdateSQ, ~1")


model.list <- read.csv(model.list.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)
model.list$model.number <- 1:nrow(model.list)
model.list

# fit the models
myobs <- ls()
myobs <- myobs[ grepl("m...",myobs,fixed=TRUE)]
cat("Removing ", myobs, "\n")
rm(list=myobs)

model.fits <- plyr::dlply(model.list, "model.number", function(x,input.data, input.ddl){
  cat("\n\n***** Starting ", unlist(x), "\n")
  #browser()
  
  #fit <- myoccMod(model=list(as.formula(paste("psi",x$psi)),
  fit <- RMark::mark(input.data, ddl=input.ddl,
                     model="Occupancy",
                     model.parameters=list(
                       Psi   =list(formula=as.formula(eval(x$Psi))),
                       p     =list(formula=as.formula(eval(x$p)))
                     )
                     #,brief=TRUE,output=FALSE, delete=TRUE
                     #,invisible=TRUE,output=TRUE  # set for debugging
  )
  mnumber <- paste("m...",formatC(x$model.number, width = 3, format = "d", flag = "0"),sep="")
  assign( mnumber, fit, envir=.GlobalEnv)
  #browser()
  fit
  
},input.data=ybf.data, input.ddl=ybf.ddl)

# examine individula model results
model.number <-4

summary(model.fits[[model.number]])
model.fits[[model.number]]$results$real
model.fits[[model.number]]$results$beta
model.fits[[model.number]]$results$derived

get.real(model.fits[[model.number]], "Psi", se=TRUE)
get.real(model.fits[[model.number]], "p",    se=TRUE)

##############################################
# Collect models
model.set <- RMark::collect.models( type="Occupancy")
model.set


names(model.set)
model.set$model.table

# Get model averaged values
Psi.ma <- RMark::model.average(model.set, param="Psi")
Psi.ma

# It is often convenient to estimate parameters (such as p) for each site using the
# values of the covariates specific for that site.
# This is similar to RPresence which gives estimates at each site.

# We create a data frame with the covariates etc as measured on the input.history dataframe.
# But we also need to add the appropriate index number for psi to the covariate data frame to account
# for the groups definition

p.covariates <- input.history
p.covariates$freq <- NULL
p.covariates$ch   <- NULL

ybf.ddl$p
p.covariates <- merge(ybf.ddl$p, p.covariates)
head(p.covariates)

# we only want a prediction for that particular model.index.
# we need to create a variable "index" that has the model.index
p.covariates$index <- p.covariates$model.index

p.pred <- covariate.predictions(model.set,
                                data=p.covariates)

# Notice that there is problem in that if model.index =1, it used jdate1 as the predictor
# and if model.index==2, it used jdate2 as the predictor. There is no easy around thi!
head(p.pred$estimates)


ggplot(data=NULL, aes( y=estimate))+
  ggtitle("Model averaged predictions of detection for each site")+
  geom_point(data=p.pred$estimates[ p.pred$estimates$index==2,], aes(x=jdate2))+
  geom_point(data=p.pred$estimates[ p.pred$estimates$index==1,], aes(x=jdate1))+
  geom_ribbon(data=p.pred$estimates[ p.pred$estimates$index==1,],aes(x=jdate1, ymin=lcl, ymax=ucl), alpha=0.2)+
  geom_ribbon(data=p.pred$estimates[ p.pred$estimates$index==2,],aes(x=jdate2, ymin=lcl, ymax=ucl), alpha=0.2)


# cleanup
cleanup(ask=FALSE)