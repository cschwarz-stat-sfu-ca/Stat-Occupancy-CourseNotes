# Multi species single season

# Co-occurence of Jordan's salamander (Plethodon jordani) (PJ) and
# members of Plethodon glutinosus (PG) in Great Smokey
# Mountains National Park (MacKenzie et al. 2004).

# 2020-06-24 CJS real$nu real$rho are now in derived$nu and derived$rho
# 2018-09-27 Code contributed by Carl James Schwarz (cschwarz.stat.sfu.cs@gmail.com)


#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------

library(ggplot2)
library(readxl)
library(reshape2)
library(RPresence)

# get the data
PG.data <- readxl::read_excel(file.path("..","Salamander2.xls"),
                              sheet="RawData", na='-',
                              col_names=FALSE,
                              range = "B3:F90")

PJ.data <- readxl::read_excel(file.path("..","Salamander2.xls"),
                              sheet="RawData", na='-',
                              col_names=FALSE,
                              range = "H3:L90")

# Convert history to compressed format
# 0=neither species detected,
# 1=only species A detected,
# 2=only species B detected,
# 3=both species detected

input.history <- PG.data + 2*PJ.data
input.history

site.covar <- readxl::read_excel(file.path("..","Salamander2.xls"),
                              sheet="RawData", col_names=TRUE,
                              range = "N2:O90")
names(site.covar)
names(site.covar) <- make.names(names(site.covar))
names(site.covar)

# Standardize Elevation covariates
elevation.mean <- mean(site.covar$Elevation..m.)
elevation.std  <- sd  (site.covar$Elevation..m.)
site.covar$Std..Elevation <- (site.covar$Elevation..m. - elevation.mean)/elevation.std



Nsites <- nrow(input.history)
Nvisits<- ncol(input.history)

# create the pao file

salamander.pao <- createPao(input.history,
                      unitcov=site.covar,
                      title="Salamander multi species - co-occurance")
summary(salamander.pao)



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


mod1 <- occMod(model=list(psi~SP+INT, p~SP+INT_o+SP:INT_o + INT_d),
               data=salamander.pao,
               type="so.2sp.1")# param="PsiBA")
summary(mod1)

names(mod1$real)

# Probability of occupancy of A
mod1$real$psiA  [1,]  

# If species occupy sites independently, then
#     psiBA = psiBa 
# i.e. regardlesss of occupancy by A, occupancy of B is 
# the same
# Pr(Occupancy of B | A present)
mod1$real$psiBA [1,]

# Pr(Occupancy of B | AA absent)
mod1$real$psiBa[1,]

# We can compute the marginal estimate of occupancy of B
psiB <-   mod1$real$psiA  * mod1$real$psiBA+
       (1-mod1$real$psiA) * mod1$real$psiBa
psiB[1,]

# look at odds ratio of occupancy of B|A to B|a.
# If there is no interaction between species then this is 1
# If less than 1, then occur less often than expected if independent
# If greater than 1, then occur more often than expected if independent
mod1$derived$nu[1,]

# For example, the marginal estimates are and prob if indendent are
c(mod1$real$psiA[1,1],psiB[1,1], mod1$real$psiA[1,1]*psiB[1,1])
# Actual estimation of occupancy of both species
mod1$real$psiBA [1,1]*mod1$real$psiA[1,1]

# Probability of detection of species A if alone
mod1$real$pA[1,]

# Probability of detection of species B if alone
mod1$real$pB[1,]

# Probability of detection when both species present
#  rA - just A; rBA B give A; rBa B given A not detected
# Again, if species do not interact, the
#      rBA = rBa
mod1$real$rA[1,]
mod1$real$rBA[1,]
mod1$real$rBa[1,]

# similarly, the odds ratio of detection of B|A vs B|a
# rho < 1, B is harder to detect if A detected compared to detecting B if A not detected
mod1$derived$rho[1,]


#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
# mod1.occind is same as model 1, but we assume that species occupancy is 
# independent. This implies that a psiBA=psiBa and is found by not
# fitting the interaction term in the occupancy mode

mod1.occind <- occMod(model=list(psi~SP,  p~SP+INT_o+INT_d+SP:INT_o),
               data=salamander.pao,
               type="so.2sp.1")
summary(mod1.occind)

names(mod1.occind$real)

# Probability of occupancy of A
mod1.occind$real$psiA  [1,]  

# If species occupy sites independently, then
#     psiBA = psiBa 
# i.e. regardlesss of occupancy by A, occupancy of B is 
# the same
# Pr(Occupancy of B | A present)
mod1.occind$real$psiBA [1,]

# Pr(Occupancy of B | AA absent)
mod1.occind$real$psiBa[1,]

# We can compute the marginal estimate of occupancy of B
psiB <-   mod1.occind$real$psiA  * mod1.occind$real$psiBA+
       (1-mod1.occind$real$psiA) * mod1.occind$real$psiBa
psiB[1,]

# look at odds ratio of occupancy of B|A to B|a.
# If there is no interaction between species then this is 1
# If less than 1, then occur less often than expected if independent
# If greater than 1, then occur more often than expected if independent
mod1.occind$derived$nu[1,]

# For example, the marginal estimates are and prob if indendent are
c(mod1.occind$real$psiA[1,1],psiB[1,1], mod1.occind$real$psiA[1,1]*psiB[1,1])
# Actual estimation of occupancy of both species
mod1.occind$real$psiBA [1,1]*mod1.occind$real$psiA[1,1]

# Probability of detection of species A if alone
mod1.occind$real$pA[1,]

# Probability of detection of species B if alone
mod1.occind$real$pB[1,]

  # Probability of detection when both species present
#  rA - just A; rBA B give A; rBa B given A not detected
# Again, if species do not interact, the
#      rBA = rBa
mod1.occind$real$rA[1,]
mod1.occind$real$rBA[1,]
mod1.occind$real$rBa[1,]

# similarly, the odds ratio of detection of B|A vs B|a
# rho < 1, B is harder to detect if A detected compared to detecting B if A not detected
mod1.occind$derived$rho[1,]


# is there any support for the indpendent occupancy model?
models<-list(mod1, 
             mod1.occind)
results<-RPresence::createAicTable(models)
summary(results)



#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
# mod1.detind is same as model 1, but we assume that detection probabilities
# of B do not depend if A has been detected or not i.e. r(B|A)=r(B|a)

mod1.detind <- occMod(model=list(psi~SP+INT,  p~SP+INT_o+SP:INT_o),
               data=salamander.pao,
               type="so.2sp.1")
summary(mod1.detind)

mod1.detind$dmat

names(mod1.detind$real)

# Probability of occupancy of A
mod1.detind$real$psiA  [1,]  

# If species occupy sites independently, then
#     psiBA = psiBa 
# i.e. regardlesss of occupancy by A, occupancy of B is 
# the same
# Pr(Occupancy of B | A present)
mod1.detind$real$psiBA [1,]

# Pr(Occupancy of B | AA absent)
mod1.detind$real$psiBa[1,]

# We can compute the marginal estimate of occupancy of B
psiB <-   mod1.detind$real$psiA  * mod1.detind$real$psiBA+
       (1-mod1.detind$real$psiA) * mod1.detind$real$psiBa
psiB[1,]

# look at odds ratio of occupancy of B|A to B|a.
# If there is no interaction between species then this is 1
# If less than 1, then occur less often than expected if independent
# If greater than 1, then occur more often than expected if independent
mod1.detind$derived$nu[1,]

# For example, the marginal estimates are and prob if indendent are
c(mod1.detind$real$psiA[1,1],psiB[1,1], mod1.detind$real$psiA[1,1]*psiB[1,1])
# Actual estimation of occupancy of both species
mod1.detind$real$psiBA [1,1]*mod1.detind$real$psiA[1,1]

# Probability of detection of species A if alone
mod1.detind$real$pA[1,]

# Probability of detection of species B if alone
mod1.detind$real$pB[1,]

# Probability of detection when both species present
#  rA - just A; rBA B give A; rBa B given A not detected
# Again, if species do not interact, the
#      rBA = rBa
mod1.detind$real$rA[1,]
mod1.detind$real$rBA[1,]
mod1.detind$real$rBa[1,]

# similarly, the odds ratio of detection of B|A vs B|a
# rho < 1, B is harder to detect if A detected compared to detecting B if A not detected
mod1.detind$derived$rho[1,]




#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
# Now do model averaging
models<-list(mod1, 
             mod1.occind,
             mod1.detind)
results<-RPresence::createAicTable(models)
summary(results)


#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
# Look at the effect of elevtion
mod.psi.elevation <- occMod(model=list(
      psi~Std..Elevation+SP+SP:Std..Elevation+INT+INT:Std..Elevation, 
      p  ~SP+INT_o+SP:INT_o),
      data=salamander.pao,
      type="so.2sp.1")
summary(mod.psi.elevation)

names(mod.psi.elevation$real)

# Probability of occupancy of A, B|A, anb B|a now depends on occupancy (on the logit scale)
# Ww need to look at the beta values
mod.psi.elevation$beta$psi
str(mod.psi.elevation$beta$psi)
mod.psi.elevation$dmat$psi

range(site.covar$Elevation..m.)
predictPSI <- data.frame(Elevation..m.=seq(400,1200,10))
predictPSI$Std..Elevation <- (predictPSI$Elevation..m. - elevation.mean)/elevation.std

predictPSI$psiA.logit  <- mod.psi.elevation$beta$psi[1,]$est + mod.psi.elevation$beta$psi[2,]$est*predictPSI$Std..Elevation
predictPSI$psiBA.logit <- mod.psi.elevation$beta$psi[3,]$est + mod.psi.elevation$beta$psi[5,]$est*predictPSI$Std..Elevation
predictPSI$psiBa.logit <- mod.psi.elevation$beta$psi[4,]$est + mod.psi.elevation$beta$psi[6,]$est*predictPSI$Std..Elevation

expit <- function (x) {1/(1+exp(-x))}
predictPSI$psiA  <- expit(predictPSI$psiA.logit)
predictPSI$psiBA <- expit(predictPSI$psiBA.logit)
predictPSI$psiBa <- expit(predictPSI$psiBa.logit)

# We can compute the marginal estimate of occupancy of B
predictPSI$psiB <-   predictPSI$psiA  * predictPSI$psiBA+
                  (1-predictPSI$psiA) * predictPSI$psiBa
predictPSI$psiAandB <- predictPSI$psiA  * predictPSI$psiBA

predictPSI[1:3,]

plotdata <- reshape2::melt(predictPSI,
                           id.vars=c("Elevation..m.","Std..Elevation"),
                           measure.vars=c("psiA","psiB","psiAandB"),
                           variable.name="Species",
                           value.name='P.Occupancy')
species.by.elev <- ggplot2::ggplot(dat=plotdata, aes(x=Elevation..m., y=P.Occupancy, color=Species))+
  ggtitle("Occupancy as a function of elevation")+
  geom_line()
species.by.elev
ggsave(plot=species.by.elev, file='species-by-elev.png',
       h=4, w=6, units="in", dpi=300)


# If species occupy sites independently, then
#     psiBA = psiBa 
# i.e. regardlesss of occupancy by A, occupancy of B is 
# the same
# Pr(Occupancy of B | A present)
mod.psi.elevation$real$psiBA [1,]

# Pr(Occupancy of B | AA absent)
mod.psi.elevation$real$psiBa[1,]

# We can compute the marginal estimate of occupancy of B
psiB <-   mod.psi.elevation$real$psiA$est  * mod.psi.elevation$real$psiBA$est+
       (1-mod.psi.elevation$real$psiA$est) * mod.psi.elevation$real$psiBa$est
psiB

# look at odds ratio of occupancy of B|A to B|a.
# If there is no interaction between species then this is 1
# If less than 1, then occur less often than expected if independent
# If greater than 1, then occur more often than expected if independent
mod.psi.elevation$derived$nu[1,]

# For example, the marginal estimates are and prob if indendent are
cbind(mod.psi.elevation$real$psiA$est,   psiB, mod.psi.elevation$real$psiA$est*psiB)
# Actual estimation of occupancy of both species
mod.psi.elevation$real$psiBA$est *mod.psi.elevation$real$psiA$est

# Probability of detection of species A if alone
mod.psi.elevation$real$pA[1,]

# Probability of detection of species B if alone
mod.psi.elevation$real$pB[1,]

# Probability of detection when both species present
#  rA - just A; rBA B give A; rBa B given A not detected
# Again, if species do not interact, the
#      rBA = rBa
mod.psi.elevation$real$rA[1,]
mod.psi.elevation$real$rBA[1,]
mod.psi.elevation$real$rBa[1,]

# similarly, the odds ratio of detection of B|A vs B|a
# rho < 1, B is harder to detect if A detected compared to detecting B if A not detected
mod.psi.elevation$derived$rho[1,]

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
# Model averaging
models<-list(mod1, 
             mod1.occind,
             mod1.detind,
             mod.psi.elevation)
results<-RPresence::createAicTable(models)
summary(results)

RPresence::modAvg(results, param="psiA")[1,]


