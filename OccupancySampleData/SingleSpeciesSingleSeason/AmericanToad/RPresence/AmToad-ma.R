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

# Fitting several models using RPresence

# 2018-08-15 Code contributed by Carl James Schwarz (cschwarz.stat.sfu.cs@gmail.com)

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

library(car)
library(readxl)
library(RPresence)
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
input.history <- input.data # the history extracted
head(input.history)

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.history)
ncol(input.history)
range(input.history, na.rm=TRUE) # check that all values are either 0 or 1
sum(is.na(input.history))    # are there any missing values?


# Get the pond information
pond.data <- readxl::read_excel("../AmericanToad.xls",
                                sheet="Pond",
                                na="-",
                                col_names=TRUE)  
head(pond.data)
pond.data$Pond <- as.factor(car::recode(pond.data$Pond,
                                        "1='P'; 0='E'; " ))
head(pond.data)


# Get the temperature data
temp.data <- readxl::read_excel("../AmericanToad.xls",
                                sheet="AmToadTemperature",
                                na="-",
                                col_names=FALSE)  # notice no column names in row 1 of data file. 
head(temp.data)
temp.data[ is.na(input.history)] <- NA

# Convert to a survey covariate. You need to stack the data by columns
survey.cov <- data.frame(Site=1:nrow(temp.data),
                         visit=as.factor(rep(1:ncol(temp.data), each=nrow(temp.data))),
                         Temperature =unlist(temp.data), stringsAsFactors=FALSE)

head(survey.cov)


# Create the *.pao file
amtoad.pao <- RPresence::createPao(input.history,
                                   unitcov=pond.data,
                                   survcov=survey.cov,
                                   title='American Toad SSSS')
amtoad.pao



# define the list of models to fit
# Notice the commas between the column and the placement of the quotes
model.list.csv <- textConnection("
p,               psi
~1,              ~1
~1,              ~Pond
~Temperature,    ~1
~Temperature,    ~Pond")

model.list <- read.csv(model.list.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)
model.list


# fit the models
model.fits <- plyr::alply(model.list, 1, function(x,detect.pao){
  cat("\n\n***** Starting ", unlist(x), "\n")
  fit <- RPresence::occMod(model=list(as.formula(paste("psi",x$psi)),
                                      as.formula(paste("p"  ,x$p  ))),
                           data=detect.pao,type="so")
  fit
},detect.pao=amtoad.pao)




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

# plot detectability as a function of temperature
head(ma.p)
ma.p$Site <- as.numeric(substring(row.names(ma.p), 4+regexpr("unit", row.names(ma.p), fixed=TRUE)))
ma.p$visit<- as.character(substr(row.names(ma.p), 2, -1+regexpr("_", row.names(ma.p), fixed=TRUE))) 

survey.cov$visit <- as.character(survey.cov$visit)

plotdata <- merge(ma.p, survey.cov)
ggplot(data=plotdata, aes(x=Temperature, y=est))+
  ggtitle("Detection probability as a function of temperature")+
  geom_point()+
  geom_ribbon(aes(ymin=lower_0.95, ymax=upper_0.95), alpha=.2)+
  ylim(0,1)
