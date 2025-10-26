# Single Species Single Season Occupancy 

## Weta Example with model averaging with a list of models to fit
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

input.data <- readxl::read_excel(file.path("..","weta.xls"),
                                    sheet="detection_histories",
                                    na="-",
                                    col_names=FALSE)  # notice no column names in row 1 of data file. 

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.data)
ncol(input.data)
range(input.data, na.rm=TRUE)
sum(is.na(input.data))

head(input.data)

# Extract the history records and create a capture history
input.history <- data.frame(freq=1,
                            ch=apply(input.data[,1:5],1,paste, collapse=""), stringsAsFactors=FALSE)
head(input.history)

# Change any NA to . in the chapter history
input.history$ch <- gsub("NA",".", input.history$ch, fixed=TRUE)
head(input.history)

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.history)


# Get the site level covariates
site_covar <- readxl::read_excel(file.path("..","weta.xls"),
                                 sheet="site_covar",
                                 na="-",
                                 col_names=TRUE)  # notice col_names in row 1 of table. 
head(site_covar)

# Create a site level covariate that is a categorical variable rather 
# than indicator variables. Be sure that it is declared as a FACTOR
input.history$BrowCat <- factor(site_covar$BrowseCat)
xtabs(~BrowCat, data=input.history,exclude=NULL, na.action=na.pass)

head(input.history)


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

Obs <-obs1*1 + obs2*2 + obs3*3
Obs<- as.data.frame(lapply(Obs, as.character), stringsAsFactors=FALSE)
Obs[ is.na(Obs)] <- "."  # set all missing values to .
head(Obs)

# make each column a factor variable. Be sure to use 
Obs[] <- lapply(Obs, as.factor)
str(Obs)

colnames(Obs) = paste("Observer", 1:5, sep="") # assign the observer at each time point
head(Obs)

Obs.df <- make.time.factor(Obs, var.name="Observer", times=1:5, intercept=1)
head(Obs.df) # Column names are indicator variables for Obsever x at time y

input.history <- cbind(input.history, Obs.df)
head(input.history)


# Illustration of using categorical covariate in the modelling process by creating groups
weta.data <- process.data(data=input.history, group="BrowCat",
                         model="Occupancy")
summary(weta.data)

# What are the parameter names for Single Season Single Species models
setup.parameters("Occupancy", check=TRUE)

# is survey covariates are present modify the ddl
weta.ddl <- make.design.data(weta.data)
weta.ddl


# define the list of models to fit
# Notice the commas between the column and the placement of the quotes
model.list.csv <- textConnection("
p,                                  psi
~1,                                 ~1
~1,                                 ~BrowCat
~BrowCat,                           ~BrowCat
~time,                              ~BrowCat
~Observer2+Observer3,               ~BrowCat
~Observer2+Observer3+time,          ~BrowCat
~Observer2+Observer3+time+BrowCat,  ~BrowCat
~Observer2+Observer3+time+BrowCat,  ~1
~Observer2+Observer3+time,          ~1")


model.list <- read.csv(model.list.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)
model.list

# fit the models
model.fits <- plyr::alply(model.list, 1, function(x,input.data, input.ddl){
  cat("\n\n***** Starting ", unlist(x), "\n")
  fit <- RMark::mark(input.data, ddl=input.ddl,
                     model="Occupancy",
                     model.parameters=list(
                       Psi   =list(formula=as.formula(eval(x$psi))),
                       p     =list(formula=as.formula(eval(x$p)))
                     ))
  fit
},input.data=weta.data, input.ddl=weta.ddl)

model.1 = model.fits[[1]]
model.2 = model.fits[[2]]
model.3 = model.fits[[3]]
model.4 = model.fits[[4]]
model.5 = model.fits[[5]]
model.6 = model.fits[[6]]
model.7 = model.fits[[7]]
model.8 = model.fits[[8]]
model.9 = model.fits[[9]]

#AICc model averaging
model.set <- collect.models(type="Occupancy")
model.set


psi.ma <- RMark::model.average(model.set, parameter="Psi", vcv=TRUE)
names(psi.ma)
head(psi.ma$estimates)


ggplot(data=psi.ma$estimates, aes(x=BrowCat, y=estimate))+
  geom_point()+
  geom_errorbar(aes(ymin=lcl, 
                    ymax=ucl), width=0.2)+
  ylim(0,1)+
  ylab("Occupancy (with 95% ci)")+xlab("Browse Category")

# or use covariate predictions to estimate the occupancy at different
# browse categories

# use the get.real to figure out the index numbers
# by looking at the alt.diff.index numbers
get.real(model.fits[[8]], parameter="Psi", se=TRUE)

psi.ma2 <- covariate.predictions(model.set, indices=11:12)
head(psi.ma2$estimates)
psi.ma2$estimates <- merge(psi.ma2$estimates, 
                           weta.ddl$Psi[,c("model.index","BrowCat")],
                           by.x="par.index", by.y="model.index")
psi.ma2$estimates

ggplot(data=psi.ma2$estimates, aes(x=BrowCat, y=estimate))+
  geom_point()+
  geom_errorbar(aes(ymin=lcl, 
                    ymax=ucl), width=0.2)+
  ylim(0,1)+
  ylab("Occupancy (with 95% ci)")+xlab("Browse Category")


p.ma <- RMark::model.average(model.set, parameter="p")
p.ma

# cleanup
cleanup(ask=FALSE)


