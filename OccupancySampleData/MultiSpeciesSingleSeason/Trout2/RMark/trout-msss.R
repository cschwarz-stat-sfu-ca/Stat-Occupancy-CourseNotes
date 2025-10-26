# Multi species single season
# Co-occurence of Bull and Brook trout
# Data provided by Parks Canada
# 2018-12-13 code submitted by Neil Faught

library(car)
library(ggplot2)
library(readxl)
library(reshape2)
library(RMark)

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

input.data <- history[,c("r1","r2")]
input.data[1:5,]

input.chistory <- data.frame(lapply(input.data, as.character), stringsAsFactors=FALSE)

input.chistory <- data.frame(lapply(input.chistory, 
                                    car::recode, " '0'='00'; '1'='10';'2'='01';'3'='11';", as.numeric=FALSE, as.factor=FALSE),
                             stringsAsFactors=FALSE)
input.history <- data.frame(ch=apply(input.chistory, 1, paste, sep="",collapse=""), freq=1)


# Change any NA to . in the chapter history
select <- grepl("NA", input.history$ch)
input.history[ select,]

input.history$ch <- gsub("NA","..", input.history$ch, fixed=TRUE)
input.history[ select,]

site.covar <- history[,c("watershed","temp","discharge")]
names(site.covar)
head(site.covar)

temp.mean <- mean(site.covar$temp)
temp.sd   <- sd  (site.covar$temp)
site.covar$temp.std <- (site.covar$temp - temp.mean)/ temp.sd

discharge.mean <- mean(site.covar$discharge)
discharge.sd   <- sd  (site.covar$discharge)
site.covar$discharge.std <- (site.covar$discharge - discharge.mean)/ discharge.sd

# Combine site covariates with input history
input.history = cbind(input.history, site.covar)
head(input.history)

# Creat RMark object
trout.data <- process.data(data=input.history,
                         model="2SpecConOccup")  # this parameterization is more stable

summary(trout.data)

# if there are time-specific covariates add them to the ddl's
trout.ddl <- RMark::make.design.data(trout.data)  # need for covariate predictions

# What are the parameter names for Single Season Single Species models
setup.parameters("2SpecConOccup", check=TRUE)

###########################################################################
# Use the psiAB parameterization # 

# Model where occupancy is not independent, and there is detection interaction
mod1 = RMark::mark(trout.data, ddl=trout.ddl,
                   model="2SpecConOccup",
                   model.parameters=list(
                     PsiA   =list(formula= ~1),
                     PsiBA  =list(formula= ~1),
                     PsiBa  =list(formula= ~1),
                     pA     =list(formula= ~1),
                     pB     =list(formula= ~1),
                     rA     =list(formula= ~1),
                     rBA    =list(formula= ~1),
                     rBa    =list(formula= ~1)
                   ))

summary(mod1)
names(mod1)

# Probability of occupancy of A
psiA = mod1$results$real[1,]
psiA

# If species occupy sites independently, then
#     psiBA = psiBa 
# i.e. regardlesss of occupancy by A, occupancy of B is 
# the same
# Pr(Occupancy of B | A present)
psiBA = mod1$results$real[2,]
psiBA

# Pr(Occupancy of B | A absent)
psiBa = mod1$results$real[3,]
psiBa

# We can compute the marginal estimate of occupancy of B
psiB = psiA[1,1]*psiBA[1,1] +
       (1-psiA[1,1])*psiBa[1,1]
psiB
# RMark can generate this too, with confindence limits
mod1$results$derived$`Occupancy of Species B`

# Look at Species Interaction Factor (SIF), which is
# similar to odds ratio of occupancy of B|A to B|a.
# SIF has the form: 
#     SIF = psi(A and B)/(psi(A)*psi(B))
#         = psi(B|A)*psi(A)/(psi(A)*(psi(A)*psi(B|A) + (1-psi(A))*psi(B|a))) 
# If there is no interaction between species then this is 1
# If less than 1, then occur less often than expected if independent
# If greater than 1, then occur more often than expected if independent
mod1$results$derived$`Species Interaction Factor`

# look at odds ratio of occupancy of B|A to B|a.
# If there is no interaction between species then this is 1
# Note, no confidence intervals (But RPrence analysis has them)
(psiBA[1,1]/(1-psiBA[1,1])) /
  (psiBa[1,1]/(1-psiBa[1,1]))

# For example, the marginal estimates are and prob if indendent are
c(psiA[1,1],psiB, psiA[1,1]*psiB)
# Actual estimation of occupancy of both species
psiBA[1,1]*psiA[1,1]

# Probability of detection of species A if alone
pA = mod1$results$real[4,]
pA

# Probability of detection of species B if alone
pB = mod1$results$real[5,]
pB

# Probability of detection when both species present
#  rA - just A; rBA B given A also detected; rBa B given A not detected
# Again, if species do not interact, then
#      rBA = rBa
rA = mod1$results$real[6,]
rA
rBA = mod1$results$real[7,]
rBA
rBa = mod1$results$real[8,]
rBa

# similarly, the odds ratio of detection of B|A vs B|a
# rho < 1, B is harder to detect if A detected compared to detecting B if A not detected
# Note: RMark doesn't automatically produce confidence interaval for rho
# RPresence does though
(rBA[1,1]/(1-rBA[1,1])) /
(rBa[1,1]/(1-rBa[1,1]))  

###########################################################################
# The detection rates are all high and almost all the same.
# Fit a simpler model where pA = pB, rBA = rBa
# Model where pA = pB = rA = rBA = rBa cannot be fit in RMark. It can be
# fit in RPresence. See http://www.phidot.org/forum/viewtopic.php?f=21&t=2974 for details
mod1.psimp = RMark::mark(trout.data, ddl=trout.ddl,
                   model="2SpecConOccup",
                   model.parameters=list(
                     PsiA   =list(formula= ~1),
                     PsiBA  =list(formula= ~1),
                     PsiBa  =list(formula= ~1),
                     pA     =list(formula= ~1, share=TRUE),
                     rA     =list(formula= ~1),
                     rBA    =list(formula= ~1, share=TRUE)
                   ))

summary(mod1.psimp)
names(mod1.psimp)

# Probability of occupancy of A
psiA = mod1.psimp$results$real[1,]
psiA

# If species occupy sites independently, then
#     psiBA = psiBa 
# i.e. regardlesss of occupancy by A, occupancy of B is 
# the same
# Pr(Occupancy of B | A present)
psiBA = mod1.psimp$results$real[2,]
psiBA

# Pr(Occupancy of B | A absent)
psiBa = mod1.psimp$results$real[3,]
psiBa

# We can compute the marginal estimate of occupancy of B
psiB = psiA[1,1]*psiBA[1,1] +
  (1-psiA[1,1])*psiBa[1,1]
psiB
# RMark can generate this too, with confindence limits
mod1.psimp$results$derived$`Occupancy of Species B`

# Look at Species Interaction Factor (SIF), which is
# similar to odds ratio of occupancy of B|A to B|a.
# SIF has the form: 
#     SIF = psi(A and B)/(psi(A)*psi(B))
#         = psi(B|A)*psi(A)/(psi(A)*(psi(A)*psi(B|A) + (1-psi(A))*psi(B|a))) 
# If there is no interaction between species then this is 1
# If less than 1, then occur less often than expected if independent
# If greater than 1, then occur more often than expected if independent
mod1.psimp$results$derived$`Species Interaction Factor`

# look at odds ratio of occupancy of B|A to B|a.
# If there is no interaction between species then this is 1
# Note, no confidence intervals (But RPrence analysis has them)
(psiBA[1,1]/(1-psiBA[1,1])) /
  (psiBa[1,1]/(1-psiBa[1,1]))

# For example, the marginal estimates are and prob if indendent are
c(psiA[1,1],psiB, psiA[1,1]*psiB)
# Actual estimation of occupancy of both species
psiBA[1,1]*psiA[1,1]

# Probability of detection of species A if alone
pA = mod1.psimp$results$real[4,]
pA

# Probability of detection of species B if alone
pB = mod1.psimp$results$real[4,]
pB

# Probability of detection when both species present
#  rA - just A; rBA B given A also detected; rBa B given A not detected
# Again, if species do not interact, then
#      rBA = rBa
rA = mod1.psimp$results$real[5,]
rA
rBA = mod1.psimp$results$real[6,]
rBA
rBa = mod1.psimp$results$real[6,]
rBa

# similarly, the odds ratio of detection of B|A vs B|a
# rho < 1, B is harder to detect if A detected compared to detecting B if A not detected
# Note: RMark doesn't automatically produce confidence interaval for rho
# RPresence does though
(rBA[1,1]/(1-rBA[1,1])) /
  (rBa[1,1]/(1-rBa[1,1]))  

# Creat AICc table for mod1 and mod1.psimp
model.set <- RMark::collect.models( type="2SpecConOccup")
model.set

rm(model.set)

###########################################################################
# mod1.occind is same as model 1, but we assume that species occupancy is 
# independent. This implies that a psiBA=psiBa.
mod1.occind = RMark::mark(trout.data, ddl=trout.ddl,
                   model="2SpecConOccup",
                   model.parameters=list(
                     PsiA   =list(formula= ~1),
                     PsiBA  =list(formula= ~1, share = TRUE),
                     pA     =list(formula= ~1),
                     pB     =list(formula= ~1),
                     rA     =list(formula= ~1),
                     rBA    =list(formula= ~1),
                     rBa    =list(formula= ~1)
                   ))

summary(mod1.occind)
names(mod1.occind)

# Probability of occupancy of A
psiA = mod1.occind$results$real[1,]
psiA

# If species occupy sites independently, then
#     psiBA = psiBa 
# i.e. regardlesss of occupancy by A, occupancy of B is 
# the same
# Pr(Occupancy of B | A present)
psiBA = mod1.occind$results$real[2,]
psiBA

# Pr(Occupancy of B | A absent)
psiBa = mod1.occind$results$real[2,]
psiBa

# We can compute the marginal estimate of occupancy of B
psiB = psiA[1,1]*psiBA[1,1] +
  (1-psiA[1,1])*psiBa[1,1]
psiB
# RMark can generate this too, with confindence limits
mod1.occind$results$derived$`Occupancy of Species B`

# Look at Species Interaction Factor (SIF), which is
# similar to odds ratio of occupancy of B|A to B|a.
# SIF has the form: 
#     SIF = psi(A and B)/(psi(A)*psi(B))
#         = psi(B|A)*psi(A)/(psi(A)*(psi(A)*psi(B|A) + (1-psi(A))*psi(B|a))) 
# If there is no interaction between species then this is 1
# If less than 1, then occur less often than expected if independent
# If greater than 1, then occur more often than expected if independent
mod1.occind$results$derived$`Species Interaction Factor`

# look at odds ratio of occupancy of B|A to B|a.
# If there is no interaction between species then this is 1
# Note, no confidence intervals (But RPrence analysis has them)
(psiBA[1,1]/(1-psiBA[1,1])) /
  (psiBa[1,1]/(1-psiBa[1,1]))

# For example, the marginal estimates are and prob if indendent are
c(psiA[1,1],psiB, psiA[1,1]*psiB)
# Actual estimation of occupancy of both species
psiBA[1,1]*psiA[1,1]

# Probability of detection of species A if alone
pA = mod1.occind$results$real[3,]
pA

# Probability of detection of species B if alone
pB = mod1.occind$results$real[4,]
pB

# Probability of detection when both species present
#  rA - just A; rBA B given A also detected; rBa B given A not detected
# Again, if species do not interact, then
#      rBA = rBa
rA = mod1.occind$results$real[5,]
rA
rBA = mod1.occind$results$real[6,]
rBA
rBa = mod1.occind$results$real[7,]
rBa

# similarly, the odds ratio of detection of B|A vs B|a
# rho < 1, B is harder to detect if A detected compared to detecting B if A not detected
# Note: RMark doesn't automatically produce confidence interaval for rho
# RPresence does though
(rBA[1,1]/(1-rBA[1,1])) /
  (rBa[1,1]/(1-rBa[1,1]))  

# Create AICc table comparing all 4 models
model.set <- RMark::collect.models( type="2SpecConOccup")
model.set

rm(model.set)

###########################################################################
# Look at the effect of temp
mod.psi.temp = RMark::mark(trout.data, ddl=trout.ddl,
                         model="2SpecConOccup",
                         model.parameters=list(
                           PsiA   =list(formula= ~temp.std),
                           PsiBA  =list(formula= ~temp.std),
                           PsiBa  =list(formula= ~temp.std),
                           pA     =list(formula= ~1, share=TRUE),
                           rA     =list(formula= ~1),
                           rBA    =list(formula= ~1, share=TRUE)
                         ))

summary(mod.psi.temp)
names(mod.psi.temp)

# Probability of occupancy of A, B|A, anb B|a now depends on occupancy (on the logit scale)
# We need to look at the beta values
mod.psi.temp$results$beta[1:6,]
str(mod.psi.temp$results$beta$estimate[1:6])

range(input.history$temp)
predictPSI <- data.frame(temp=seq(3,12,.1))
predictPSI$temp.std <- (predictPSI$temp - temp.mean)/temp.sd

predictPSI$psiA.logit  <- mod.psi.temp$results$beta$estimate[1] + 
                          mod.psi.temp$results$beta$estimate[2]*predictPSI$temp.std
predictPSI$psiBA.logit <- mod.psi.temp$results$beta$estimate[3] + 
                          mod.psi.temp$results$beta$estimate[4]*predictPSI$temp.std
predictPSI$psiBa.logit <- mod.psi.temp$results$beta$estimate[5] + 
                          mod.psi.temp$results$beta$estimate[6]*predictPSI$temp.std

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


# Probability of occupancy of A
psiA = mod.psi.temp$results$real[1,]
psiA

# If species occupy sites independently, then
#     psiBA = psiBa 
# i.e. regardlesss of occupancy by A, occupancy of B is 
# the same
# Pr(Occupancy of B | A present)
psiBA = mod.psi.temp$results$real[2,]
psiBA

# Pr(Occupancy of B | A absent)
psiBa = mod.psi.temp$results$real[3,]
psiBa

# We can compute the marginal estimate of occupancy of B
psiB = psiA[1,1]*psiBA[1,1] +
  (1-psiA[1,1])*psiBa[1,1]
psiB
# RMark can generate this too, with confindence limits
mod.psi.temp$results$derived$`Occupancy of Species B`

# Look at Species Interaction Factor (SIF), which is
# similar to odds ratio of occupancy of B|A to B|a.
# SIF has the form: 
#     SIF = psi(A and B)/(psi(A)*psi(B))
#         = psi(B|A)*psi(A)/(psi(A)*(psi(A)*psi(B|A) + (1-psi(A))*psi(B|a))) 
# If there is no interaction between species then this is 1
# If less than 1, then occur less often than expected if independent
# If greater than 1, then occur more often than expected if independent
mod.psi.temp$results$derived$`Species Interaction Factor`

# look at odds ratio of occupancy of B|A to B|a.
# If there is no interaction between species then this is 1
# Note, no confidence intervals (But RPrence analysis has them)
(psiBA[1,1]/(1-psiBA[1,1])) /
  (psiBa[1,1]/(1-psiBa[1,1]))

# For example, the marginal estimates are and prob if indendent are
c(psiA[1,1],psiB, psiA[1,1]*psiB)
# Actual estimation of occupancy of both species
psiBA[1,1]*psiA[1,1]

# Probability of detection of species A if alone
pA = mod.psi.temp$results$real[4,]
pA

# Probability of detection of species B if alone
pB = mod.psi.temp$results$real[4,]
pB

# Probability of detection when both species present
#  rA - just A; rBA B given A also detected; rBa B given A not detected
# Again, if species do not interact, then
#      rBA = rBa
rA = mod.psi.temp$results$real[5,]
rA
rBA = mod.psi.temp$results$real[6,]
rBA
rBa = mod.psi.temp$results$real[6,]
rBa

# similarly, the odds ratio of detection of B|A vs B|a
# rho < 1, B is harder to detect if A detected compared to detecting B if A not detected
# Note: RMark doesn't automatically produce confidence interaval for rho
# RPresence does though
(rBA[1,1]/(1-rBA[1,1])) /
  (rBa[1,1]/(1-rBa[1,1]))  


###########################################################################


# Create AICc table comparing all 4 models
model.set <- RMark::collect.models( type="2SpecConOccup")
model.set

# No need for model averaging as 1 model contains 100% of weight

# cleanup
cleanup(ask=FALSE)
