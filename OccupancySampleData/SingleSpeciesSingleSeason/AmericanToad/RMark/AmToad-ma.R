# American Toad

# Extracted from
#    Darryl I. MacKenzie, et al. 2002.
#    Estimating site occupancy rates when detection probabilities are less than one.
#    Ecology 83:2248-2255.
#    doi:10.1890/0012-9658(2002)083[2248:ESORWD]2.0.CO;2]

# 29 sites with 82 sampling occasions in 2000.
# Volunteers visited sites and recorded presence/absence of toads by calls.
# Habitat (type of pond, permanent or ephemeral) and temperature at visit recorded.

# Single Species Single Season Occupancy

# Fitting models using RMark

# 2018-08-15 Code contributed by Carl James Schwarz (cschwarz.stat.sfu.cs@gmail.com)


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

input.data <- readxl::read_excel("../AmericanToad.xls",
                                 sheet="AmToadDetectionHistories",
                                 na="-",
                                 col_names=FALSE)  # notice no column names in row 1 of data file. 

head(input.data)


# Extract the history records

# Extract the history records and create a capture history
input.history <- data.frame(freq=1,
                            ch=apply(input.data,1,paste, collapse=""), stringsAsFactors=FALSE)
head(input.history)

# Change any NA to . in the chapter history
input.history$ch <- gsub("NA",".", input.history$ch, fixed=TRUE)
head(input.history)

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.history)




# Get the pond information
pond.data <- readxl::read_excel("../AmericanToad.xls",
                                sheet="Pond",
                                na="-",
                                col_names=TRUE)  
head(pond.data)
input.history$Pond <- as.factor(car::recode(pond.data$Pond,
                                        "1='P'; 0='E'; " ))
head(input.history)




# Get the temperature data
temp.data <- readxl::read_excel("../AmericanToad.xls",
                                sheet="AmToadTemperature",
                                na="-",
                                col_names=FALSE)  # notice no column names in row 1 of data file. 
head(temp.data)
# Notice that RMark does not allow missing values in time-varying covariates, even when visits are not made
# so do not set these to missing as in RPresence
# temp.data[ is.na(input.data)] <- NA


colnames(temp.data) = paste("Temp", 1:ncol(temp.data), sep="") # assign the observer at each time point
head(temp.data)

input.history <- cbind(input.history, temp.data)
head(input.history)


# Illustration of using categorical covariate in the modelling process by creating groups
amtoad.data <- process.data(data=input.history, group="Pond",
                         model="Occupancy")
summary(amtoad.data)

# Is survey covariates are to be used for p, modify the ddl
amtoad.ddl <- make.design.data(amtoad.data)
amtoad.ddl

# What are the parameter names for Single Season Single Species models
setup.parameters("Occupancy", check=TRUE)



# define the list of models to fit
# Notice the commas between the column and the placement of the quotes
model.list.csv <- textConnection("
p,         Psi
~1,        ~1
~1,        ~Pond
~Temp,    ~1
~Temp,    ~Pond")


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
  
},input.data=amtoad.data, input.ddl=amtoad.ddl)


# examine individula model results
model.number <-2

summary(model.fits[[model.number]])
model.fits[[model.number]]$results$real
model.fits[[model.number]]$results$beta
model.fits[[model.number]]$results$derived

get.real(model.fits[[model.number]], "Psi", se=TRUE)
get.real(model.fits[[model.number]], "p",    se=TRUE)


# collect models and make AICc table

model.set <- RMark::collect.models( type="Occupancy")
model.set

# cleanup
cleanup()


