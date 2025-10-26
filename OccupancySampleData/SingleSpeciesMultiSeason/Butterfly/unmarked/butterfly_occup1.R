rm(list=ls())


#install unmarked if not installed
#install.packages("unmarked")
#load library
library(unmarked)
bfly<-read.csv("../Butterfly.csv")
bfly<-subset(bfly,select=c(Field,Border,Transect,Treatment,Date,Visit,CLLESSU))
library(reshape)
bfly$Transect=as.factor(bfly$Transect)

junk<-melt(bfly,id.var=c("Transect","Visit"),measure.var="CLLESSU")

j2=cast(junk,Transect ~ Visit)

covs<-subset(bfly,select=c(Transect,Treatment))
covs<-as.data.frame(with(covs,table(Transect,Treatment)))
covs<-subset(covs,Freq>0)

site.cov<-merge(j2,covs,by="Transect")[,8]
site.cov<-data.frame(treatment=site.cov)


y<-j2[,2:7]

y<-cbind(j2$"1",j2$"2",j2$"3",j2$"4",j2$"5",j2$"6")

S <- nrow(y) # number of sites
J <- 2 # number of secondary sampling occasions
T <- 3 # number of primary period
#dummy variables to create seasonal effects (need 2 for 3 occasions)
season<-as.factor(rep(c(1,2,3),S))

yearly.cov<-data.frame(season)
#additive time effects within season. to model effects change between seasons create an interaction term
secondary<-as.factor(rep(c(1,2,1,2,1,2),S))

second.cov<-data.frame(secondary)
butterfly<-unmarkedMultFrame(y=y,siteCovs=site.cov, numPrimary=T,yearlySiteCovs=yearly.cov,obsCovs=second.cov)
fm1<-colext(psiformula = ~1, gammaformula = ~1, epsilonformula = ~1, 
    pformula = ~1, data=butterfly) 
#back transformations
backTransform(fm1,'det')
backTransform(fm1,"psi")
backTransform(fm1,"col")
backTransform(fm1,"ext")
fm2<-colext(psiformula = ~treatment, gammaformula = ~1, epsilonformula = ~1, 
    pformula = ~1, data=butterfly) 
fm3<-colext(psiformula= ~1, gammaformula = ~treatment, epsilonformula = ~1, 
    pformula = ~1, data=butterfly) 

fm4<-colext(psiformula= ~1, gammaformula = ~1, epsilonformula = ~treatment, 
    pformula = ~1, data=butterfly) 

fm1
fm2
fm3
fm4