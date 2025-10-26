# Multi species single season

# Co-occurence of Bull and Brook trout
# Data provided by Parks Canada

# 2020-06-24 CJS real$rho and real$nu now in derived$

#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------

library(ggplot2)
library(plyr)
library(readxl)
library(reshape2)
library(RPresence)

# get the data
occ.data <- readxl::read_excel("../Trout2.xlsx",
                               sheet="fish_sp_stacked_juv_delta")
head(occ.data)

# check site code
xtabs(~site+rep, data=occ.data, exclude=NULL, na.action=na.pass)

# some site numbers are "reused", but a combination of watershed x site is unique
nmeasure <- plyr::ddply(occ.data, c("watershed","site"), plyr::summarize, count=length(rep))
nmeasure[ nmeasure$count != 4,]

xtabs(~presence, data=occ.data, exclude=NULL, na.action=na.pass)
occ.data$presence <- as.numeric(occ.data$presence)

# create the detection records for each rep
# Species A = BKTR; Species B=BLTR
# it is believed that BKTR have a preference for warmer water (see later)
# so we code them as the "dominant" species (Species A) so that 
# a linear logisitic regression makes sense for it.
# If you change the coding, the program doesn't coverge well (perhaps a 
# change in initial values may be helpful)
dhist1 <- plyr::ddply(occ.data, c("watershed","site","rep"), plyr::summarize,
                      history=sum(presence*1*(species=="BKTR")+
                                  presence*2*(species=="BLTR")),
                      temp = mean(temp),
                      discharge=mean(discharge))

history <- reshape2::dcast(dhist1, watershed+site+temp+discharge~rep,
                    variable.name="rep",
                    value.var="history")
head(history)
history <- plyr::rename(history, c("1"="r1", "2"="r2"))
xtabs(~r1+r2, data=history, exclude=NULL, na.action=na.pass)

input.history <- history[,c("r1","r2")]
input.history[1:5,]

site.covar <- history[,c("watershed","temp","discharge")]
names(site.covar)
head(site.covar)

temp.mean <- mean(site.covar$temp)
temp.sd   <- sd  (site.covar$temp)
site.covar$temp.std <- (site.covar$temp - temp.mean)/ temp.sd

discharge.mean <- mean(site.covar$discharge)
discharge.sd   <- sd  (site.covar$discharge)
site.covar$discharge.std <- (site.covar$discharge - discharge.mean)/ discharge.sd


Nsites <- nrow(input.history)
Nvisits<- ncol(input.history)

# create the pao file

trout.pao <- createPao(input.history,
                      unitcov=site.covar,
                      title="trout multi species - co-occurance")
summary(trout.pao)

#-----
# Use the psiAB parameterization 

mod1 <- occMod(model=list(psi~SP+INT, p~SP+INT_o+SP:INT_o + INT_d),
               data=trout.pao,
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

# Pr(Occupancy of B | A absent)
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
#  rA - just A; rBA B given A also detected; rBa B given A not detected
# Again, if species do not interact, then
#      rBA = rBa
mod1$real$rA[1,]
mod1$real$rBA[1,]
mod1$real$rBa[1,]

# similarly, the odds ratio of detection of B|A vs B|a
# rho < 1, B is harder to detect if A detected compared to detecting B if A not detected
mod1$derived$rho[1,]


#------
# The detection rates are all high and almost all the same.
# Fit a much simpler model

mod1.psimp <- occMod(model=list(psi~SP+INT, p~1),
               data=trout.pao,
               type="so.2sp.1")# param="PsiBA")
summary(mod1.psimp)

names(mod1.psimp$real)

# Probability of occupancy of A
mod1.psimp$real$psiA  [1,]  

# If species occupy sites independently, then
#     psiBA = psiBa 
# i.e. regardlesss of occupancy by A, occupancy of B is 
# the same
# Pr(Occupancy of B | A present)
mod1.psimp$real$psiBA [1,]

# Pr(Occupancy of B | AA absent)
mod1.psimp$real$psiBa[1,]

# We can compute the marginal estimate of occupancy of B
psiB <-   mod1.psimp$real$psiA  * mod1.psimp$real$psiBA+
       (1-mod1.psimp$real$psiA) * mod1.psimp$real$psiBa
psiB[1,]

# look at odds ratio of occupancy of B|A to B|a.
# If there is no interaction between species then this is 1
# If less than 1, then occur less often than expected if independent
# If greater than 1, then occur more often than expected if independent
mod1.psimp$derived$nu[1,]

# For example, the marginal estimates are and prob if indendent are
c(mod1.psimp$real$psiA[1,1],psiB[1,1], mod1.psimp$real$psiA[1,1]*psiB[1,1])
# Actual estimation of occupancy of both species
mod1.psimp$real$psiBA [1,1]*mod1.psimp$real$psiA[1,1]

# Probability of detection of species A if alone
mod1.psimp$real$pA[1,]

# Probability of detection of species B if alone
mod1.psimp$real$pB[1,]

# Probability of detection when both species present
#  rA - just A; rBA B give A; rBa B given A not detected
# Again, if species do not interact, the
#      rBA = rBa
mod1.psimp$real$rA[1,]
mod1.psimp$real$rBA[1,]
mod1.psimp$real$rBa[1,]

# similarly, the odds ratio of detection of B|A vs B|a
# rho < 1, B is harder to detect if A detected compared to detecting B if A not detected
mod1.psimp$derived$rho[1,]

models<-list(mod1, 
             mod1.psimp)
results<-RPresence::createAicTable(models)
summary(results)


#------------------
# mod1.occind is same as model 1, but we assume that species occupancy is 
# independent. This implies that a psiBA=psiBa and is found by not
# fitting the interaction term in the occupancy mode

mod1.occind <- occMod(model=list(psi~SP,  p~1),
               data=trout.pao,
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




models<-list(mod1,
             mod1.psimp,
             mod1.occind)
results<-RPresence::createAicTable(models)
summary(results)


#-----
# Look at the effect of temp
mod.psi.temp <- occMod(model=list(
      psi~temp.std+SP+SP:temp.std+INT+INT:temp.std, 
      p  ~1),
      data=trout.pao,
      type="so.2sp.1")
summary(mod.psi.temp)

names(mod.psi.temp$real)

# Probability of occupancy of A, B|A, anb B|a now depends on occupancy (on the logit scale)
# We need to look at the beta values
mod.psi.temp$beta$psi
str(mod.psi.temp$beta$psi)
mod.psi.temp$dmat$psi

range(site.covar$temp)
predictPSI <- data.frame(temp=seq(3,12,.1))
predictPSI$temp.std <- (predictPSI$temp - temp.mean)/temp.sd

predictPSI$psiA.logit  <- mod.psi.temp$beta$psi[1,1] + mod.psi.temp$beta$psi[2,1]*predictPSI$temp.std
predictPSI$psiBA.logit <- mod.psi.temp$beta$psi[3,1] + mod.psi.temp$beta$psi[5,1]*predictPSI$temp.std
predictPSI$psiBa.logit <- mod.psi.temp$beta$psi[4,1] + mod.psi.temp$beta$psi[6,1]*predictPSI$temp.std

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
                           id.vars=c("temp","temp.std"),
                           measure.vars=c("psiA","psiB","psiAandB"),
                           variable.name="Species",
                           value.name='P.Occupancy')
species.by.temp <- ggplot2::ggplot(dat=plotdata, aes(x=temp, y=P.Occupancy, color=Species))+
  ggtitle("Occupancy as a function of temp\nSpecies A = BKTR; Species B=BLTR")+
  geom_line()
species.by.temp
ggsave(plot=species.by.temp, file='species-by-temp.png',
       h=4, w=6, units="in", dpi=300)


# If species occupy sites independently, then
#     psiBA = psiBa 
# i.e. regardlesss of occupancy by A, occupancy of B is 
# the same
# Pr(Occupancy of B | A present)
mod.psi.temp$real$psiBA [1,]

# Pr(Occupancy of B | A absent)
mod.psi.temp$real$psiBa[1,]

# We can compute the marginal estimate of occupancy of B
psiB <-   mod.psi.temp$real$psiA  * mod.psi.temp$real$psiBA+
       (1-mod.psi.temp$real$psiA) * mod.psi.temp$real$psiBa
psiB[1,]

# look at odds ratio of occupancy of B|A to B|a.
# If there is no interaction between species then this is 1
# If less than 1, then occur less often than expected if independent
# If greater than 1, then occur more often than expected if independent
mod.psi.temp$derived$nu[1,]

# For example, the marginal estimates are and prob if indendent are
c(mod.psi.temp$real$psiA[1,1],psiB[1,1], mod.psi.temp$real$psiA[1,1]*psiB[1,1])
# Actual estimation of occupancy of both species
mod.psi.temp$real$psiBA [1,1]*mod.psi.temp$real$psiA[1,1]

# Probability of detection of species A if alone
mod.psi.temp$real$pA[1,]

# Probability of detection of species B if alone
mod.psi.temp$real$pB[1,]

# Probability of detection when both species present
#  rA - just A; rBA B give A; rBa B given A not detected
# Again, if species do not interact, the
#      rBA = rBa
mod.psi.temp$real$rA[1,]
mod.psi.temp$real$rBA[1,]
mod.psi.temp$real$rBa[1,]

# similarly, the odds ratio of detection of B|A vs B|a
# rho < 1, B is harder to detect if A detected compared to detecting B if A not detected
mod.psi.temp$derived$rho[1,]


#-------------
# Model averaging
models<-list(mod1, 
             mod1.occind,
             mod1.psimp,
             mod.psi.temp)
results<-RPresence::createAicTable(models)
summary(results)

RPresence::modAvg(results, param="psiA")[1,]


