# Multi species single season with a list of models

# Using RPresence

# Co-occurence of spotted (SO) and barred owls (BO)

# 2020-06-24 CJS real$rho and real$nu now in derived$
# 2018-12-13 code contributed by Neil Faught

#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------

library(car)
library(ggplot2)
library(readxl)
library(reshape2)
library(RPresence)

# get the data
SO.data <- readxl::read_excel(file.path("..","SpottedAndBarredOwls.xls"),
                              sheet="RawData", na='-',
                              col_names=FALSE,
                              range = "B11:K161")

BO.data <- readxl::read_excel(file.path("..","SpottedAndBarredOwls.xls"),
                              sheet="RawData", na='-',
                              col_names=FALSE,
                              range = "M11:V161")

# Convert history to compressed format
# 0=neither species detected,
# 1=only species A detected,
# 2=only species B detected,
# 3=both species detected

input.data <- 2*BO.data + SO.data
input.data

input.history = input.data

# Read in survey level covariates
night = readxl::read_excel(file.path("..","SpottedAndBarredOwls.xls"),
                           sheet="Nite Covariate", na='-',
                           col_names=FALSE,
                           range = "B11:K161")
head(night)

# Convert to a survey covariate. You need to stack the data by columns
survey.cov <- data.frame(Site=1:nrow(night),
                         visit=as.factor(rep(1:ncol(night), each=nrow(night))),
                         night =unlist(night), stringsAsFactors=FALSE)

head(survey.cov)

# create the pao file

owl.pao <- createPao(input.history,
                            survcov = survey.cov,
                            title="Owl multi species - co-occurance")
summary(owl.pao)

#-----
# Use the psiAB parameterization 
#    Parameters of the model are 
#        psiA  - occupancy probability of species A
#        psiBA - occupancy probability of species B if species A is present
#        psiBa - occupancy probability of species B if species A is absent
#
#    If species are independent thatn psiBA = psiBa.
#       Alternatively, nu = odds(B|A)/odds(B|a) = 1.
#
#    Detection parameters
#        pA    - probability of detection if A is alone in the site
#        pB    - probability of detection if B is alone in the site
#        rA    - probability of detecting A only given both on the site
#        rBA   - probability of detecting B given that A was detected and both are on the site
#        rBa   - probability of detecting B given that A not detected and both are on the site 
#    Ifspecies do not interact, then
#        rBA = rBa
#    and
#        rho = odds(rBA)/odds(rBa) = 1


#
# The following special variables are available for modelling  PSI
#    SP species effect 
#    INT interaction of effect on presence of species B when species A was, or was not present 
#
# Model for PSI.... impact on parameters
#    psi~1      psiA=psiBA=psiBa       (1 parameter)
#    psi~SP     psiA    psiBA=psiBa    (2 parameters)
#    psi~SP+INT psiA    psiBA   psiBa  (3 parameters)
#
# The following special variables are available for p
#   SP species effect
#   INT_o is a detection-level interaction where the OCCUPANCY of one species 
#         changes the detection probability of the other species 
#   INT_d is a detection-level interaction where DETECTION of one species changes the 
#         detection probability of the other species in the same survey. 
#
# Model for p.... impact on parameters
#    p~1                          pA = pB = rA = rBA = rBa   (1 parameter)
#    p~SP                         pA=rA,  pB=rBA=rBa         (2 parameters)
#    p~SP+INT_o                   pA=rA,  pB, rBA=rBa        (3 parameters)
#    p~SP+INT_o+SP:INT_o          pA, rA, pB, rBA=rBa        (4 parameters) 
#    p~SP+INT_o+SP:INT_o+INT_d    pA, rA, pB, rBA, rBa       (5 parameters) 




# define the list of models to fit
# Notice the commas between the column and the placement of the quotes
model.list.csv <- textConnection("
p,                               psi
~SP+INT_o+SP:INT_o + INT_d,              ~SP+INT
~SP+INT_o+SP:INT_o,                      ~SP+INT
~SP+INT_o+SP:INT_o+ INT_d,               ~SP
~SP+INT_o+SP:INT_o,                      ~SP
~SP+INT_o + SP:INT_o + SP:INT_o + night, ~SP+INT")

model.list <- read.csv(model.list.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)
model.list


# fit the models
model.fits <- plyr::alply(model.list, 1, function(x,detect.pao){
  cat("\n\n***** Starting ", unlist(x), "\n")
  fit <- RPresence::occMod(model=list(as.formula(paste("psi",x$psi)),
                                      as.formula(paste("p"  ,x$p  ))),
                           data=detect.pao,type="so.2sp.1", randint=10)
  fit
},detect.pao=owl.pao)



# Look the output from a specific model
check.model <- 5

names(model.fits[[check.model]])
model.fits[[check.model]]$beta
model.fits[[check.model]]$dmat

names(model.fits[[check.model]]$real)
model.fits[[check.model]]$real$psiA[1:5,]
model.fits[[check.model]]$real$psiBA[1:5,]
model.fits[[check.model]]$real$psiBa[1:5,]
model.fits[[check.model]]$derived$nu[1:5,]

model.fits[[check.model]]$real$pA[1:5,]
model.fits[[check.model]]$real$pB[1:5,]
model.fits[[check.model]]$real$rA[1:5,]
model.fits[[check.model]]$real$rBA[1:5,]
model.fits[[check.model]]$real$rBa[1:5,]
model.fits[[check.model]]$derived$rho[1:5,]

names(model.fits[[check.model]]$derived)




# Model averaging
aic.table <- RPresence::createAicTable(model.fits)
aic.table$table
