# Multi species single season with a list of models

# Co-occurence of Grizzly and black bears as measured by camera
#
# Data from
#   Ladle A, Steenweg R, Shepherd B, Boyce MS (2018)
#   The role of human outdoor recreation in shaping patterns of grizzly bear-black 
#   bear co-occurrence. 
#   PLOS ONE 13(2): e0191730. https://doi.org/10.1371/journal.pone.0191730

#   Ladle A, Steenweg R, Shepherd B, Boyce MS (2018)
#   Data from: The role of human outdoor recreation in shaping patterns of 
#   grizzly bear-black bear co-occurrence. 
#   Dryad Digital Repository. https://doi.org/10.5061/dryad.s81k3

# Briefly,
# Cameras placed on trails in summer 2014
# Classified images as type of bear and type of human activity.

# 2020-06-24 CJS real$rho and real$nu now in derived$
# 2018-12-26 Code contributed by Carl James Schwarz (cschwarz.stat.sfu.cs@gmail.com)


#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------

library(ggplot2)
library(readxl)
library(reshape2)
library(RPresence)

# get the data. Sites x # of detections grouped into 4 day intervals
gb.data <- read.csv(file.path("..","Raw data","four_day_grizzly.csv"), header=FALSE, strip.white=TRUE, as.is=TRUE)
bb.data <- read.csv(file.path("..","Raw data","four_day_black.csv"),   header=FALSE, strip.white=TRUE, as.is=TRUE)

# Notice that in some days, more than 1 detection recorded.
range(gb.data, na.rm=TRUE)
apply(gb.data, 1, max, na.rm=TRUE)

# Convert history to compressed format
# 0=neither species detected,
# 1=only species A detected,
# 2=only species B detected,
# 3=both species detected

# we need to ignore multple detections of each species in a camera
input.history <- 1*(bb.data>=1) + 2*(gb.data>=1)
input.history


# Get the site covariates

site.covar <- read.csv(file.path("..","Raw data","covariates.csv"),   header=TRUE, strip.white=TRUE, as.is=TRUE)
names(site.covar)
names(site.covar) <- make.names(names(site.covar))
names(site.covar)

head(site.covar)


# Standardize Elevation covariates
elev.mean <- mean(site.covar$elev)
elev.s    <- sd  (site.covar$elev)
site.covar$elev.std <- (site.covar$elev - elev.mean)/elev.s



Nsites <- nrow(input.history)
Nvisits<- ncol(input.history)

# create the pao file

bear.pao <- createPao(input.history,
                            unitcov=site.covar,
                            title="bear multi species - co-occurance")
summary(bear.pao)



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
~SP+INT_o+SP:INT_o + INT_d,              ~SP
~SP+INT_o+SP:INT_o  ,                    ~SP+INT
~SP+INT_o+SP:INT_o,                      ~SP + elev.std+SP:elev.std+INT+INT:elev.std")

model.list <- read.csv(model.list.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)
model.list


# fit the models
model.fits <- plyr::alply(model.list, 1, function(x,detect.pao){
  cat("\n\n***** Starting ", unlist(x), "\n")
  fit <- RPresence::occMod(model=list(as.formula(paste("psi",x$psi)),
                                      as.formula(paste("p"  ,x$p  ))),
                           data=detect.pao,type="so.2sp.1", randint=10)
  fit
},detect.pao=bear.pao)




# Look the output from a specific model
check.model <- 4

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

names(aic.table)



ma.psiA<- RPresence::modAvg(aic.table, param="psiA")
ma.psiA$parameter <- "psiA"
ma.psiA <- cbind(ma.psiA, site.covar) 
  
ma.psiBA<- RPresence::modAvg(aic.table, param="psiBA")
ma.psiBA$parameter <- "psiBA"
ma.psiBA <- cbind(ma.psiBA, site.covar) 

ma.psiBa<- RPresence::modAvg(aic.table, param="psiBa")
ma.psiBa$parameter <- "psiBa"
ma.psiBa <- cbind(ma.psiBa, site.covar) 

ma.psiB <- site.covar
ma.psiB$parameter <- 'psiB'
ma.psiB$est <- ma.psiBA$est * ma.psiA$est +
               ma.psiBa$est * (1-ma.psiA$est)
# too hard to compute the se

psi.all <- plyr::rbind.fill(ma.psiA, ma.psiBA, ma.psiBa, ma.psiB)

head(psi.all)

ggplot(data=psi.all, aes(x=elev, y=est))+
   ggtitle("Estimated occupancy parameters vs. elevation")+
   geom_point()+
   geom_ribbon( aes(ymin=lower_0.95, ymax=upper_0.95),alpha=0.2)+
   facet_wrap(~parameter, ncol=2)

# Why is psiBa so bad? This may be a case of complete separation.
# The detection probabilities are all fairly large so that over 5 visits, there is almost certain detection
# if the species is present. Consequently, the observed occupancy is likely a pretty good indication of the 
# actual occupancy.
RPresence::modAvg(aic.table, param="pA")[1,]
RPresence::modAvg(aic.table, param="pB")[1,]
RPresence::modAvg(aic.table, param="rA")[1,]
RPresence::modAvg(aic.table, param="rBA")[1,]
RPresence::modAvg(aic.table, param="rBa")[1,]

# Compute the observed occupance of each species at each camera
# Don't forget that cameras record number of times bear are seen
obs.occ <- data.frame(elev = site.covar$elev,
                      oogb = pmin(1, apply(gb.data,1,max, na.rm=TRUE)),
                      oobb = pmin(1, apply(bb.data,1,max, na.rm=TRUE)))
head(obs.occ)

ggplot(data=obs.occ, aes(x=elev, y=oobb))+
   ggtitle("Observed occupancy of bb given gb")+
   geom_point()+
   geom_smooth(method="glm", method.args = list(family = quasibinomial(link = 'logit')))+
   facet_wrap(~oogb, ncol=1)

# and look at when the species are reversed
ggplot(data=obs.occ, aes(x=elev, y=oogb))+
   ggtitle("Observed occupancy of gb given bb")+
   geom_point()+
   geom_smooth(method="glm", method.args = list(family = quasibinomial(link = 'logit')))+
   facet_wrap(~oobb, ncol=1)






