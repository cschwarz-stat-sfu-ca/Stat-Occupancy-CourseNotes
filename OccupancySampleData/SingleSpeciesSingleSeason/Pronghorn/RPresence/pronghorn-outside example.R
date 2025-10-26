# load required R libraries
library(RPresence) # not available on CRAN. https://www.mbr-pwrc.usgs.gov/software/presence.shtml
library(ggplot2)
library(AICcmodavg)


# load in data file
# open this file in excel to examine it's structure if you're not familiar with R  
pronghorn<-read.csv(file.path("..","pronghorn.csv"), header=TRUE, as.is=TRUE, strip.white=TRUE)

# slope covariate in degrees. Convert to proportion slope
pronghorn$slope<-pronghorn$slope/90

# distance to water in m. Convert to km
pronghorn$DW<-pronghorn$DW/1000

## create pao data object. Survey.1 and Survey.2 is the detection/nondetection data
prong<-createPao(data=pronghorn[,c("Survey.1","Survey.2")],
                 unitcov = pronghorn[,4:7]) 

###############################
# prepare some data frames that will be used later for creating graphs 
sagebrush<-expand.grid(sagebrush=seq(0,20,length.out=100),slope=median(prong$unitcov$slope),DW=median(prong$unitcov$DW),aspect=c("N","E","S","W"))
slope<-expand.grid(sagebrush=median(prong$unitcov$sagebrush),slope=seq(0,1,length.out=100),DW=median(prong$unitcov$DW),aspect=c("N","E","S","W"))
DW<-expand.grid(sagebrush=median(prong$unitcov$sagebrush),slope=median(prong$unitcov$slope),DW=seq(0,3,length.out=100),aspect=c("N","E","S","W"))
aspect<-expand.grid(sagebrush=median(prong$unitcov$sagebrush),slope=median(prong$unitcov$slope),DW=median(prong$unitcov$DW),aspect=c("N","E","S","W"))

#################################################################################################
### logistic regression using GLM
### 'naive' analysis that doesn't account for detection probability
#################################################################################################
glm.data<-prong$unitcov

# y =1 if ever detected at plot, =0 if never detected
glm.data$y<-as.numeric(rowSums(prong$det)>0) 

glm.mod1<-glm(y~1,data=glm.data,family=binomial)
glm.mod2<-glm(y~sagebrush,data=glm.data,family=binomial)
glm.mod3<-glm(y~slope,data=glm.data,family=binomial)
glm.mod4<-glm(y~sagebrush+slope,data=glm.data,family=binomial)
glm.mod5<-glm(y~DW,data=glm.data,family=binomial)
glm.mod6<-glm(y~sagebrush+DW,data=glm.data,family=binomial)
glm.mod7<-glm(y~slope+DW,data=glm.data,family=binomial)
glm.mod8<-glm(y~sagebrush+slope+DW,data=glm.data,family=binomial)
glm.mod9<-glm(y~aspect,data=glm.data,family=binomial)
glm.mod10<-glm(y~sagebrush+aspect,data=glm.data,family=binomial)
glm.mod11<-glm(y~slope+aspect,data=glm.data,family=binomial)
glm.mod12<-glm(y~sagebrush+slope+aspect,data=glm.data,family=binomial)
glm.mod13<-glm(y~DW+aspect,data=glm.data,family=binomial)
glm.mod14<-glm(y~sagebrush+DW+aspect,data=glm.data,family=binomial)
glm.mod15<-glm(y~slope+DW+aspect,data=glm.data,family=binomial)
glm.mod16<-glm(y~sagebrush+slope+DW+aspect,data=glm.data,family=binomial)

mods<-list(glm.mod1,glm.mod2,glm.mod3,glm.mod4,
           glm.mod5,glm.mod6,glm.mod7,glm.mod8,
           glm.mod9,glm.mod10,glm.mod11,glm.mod12,
           glm.mod13,glm.mod14,glm.mod15,glm.mod16)

## AIC comparison of GLM models
aictab(mods,second.ord=FALSE)
summary(mods[[7]]) # top aic-ranked model

#####################################################
# produce some plots of model-averaged estimates to
# illustrate estimated effects
#####################################################

# predict presence probability for different sagebrush values, while fixing other values
glm.ma.sagebrush <- modavgPred(mods,type="response",newdata=sagebrush,uncond.se="revised")

## manually calculate confidence intervals
logit<-qlogis(glm.ma.sagebrush$mod.avg.pred)
l.se<-glm.ma.sagebrush$uncond.se/(glm.ma.sagebrush$mod.avg.pred*(1-glm.ma.sagebrush$mod.avg.pred))
glm.ma.sagebrush$lower<-plogis(logit-1.96*l.se)
glm.ma.sagebrush$upper<-plogis(logit+1.96*l.se)

#jpeg("PronghornNaiveSagebrush.jpg",width=800,height=800,res=144)
#Convert model averaged results to a data frame
glm.ma.sagebrush = as.data.frame(glm.ma.sagebrush)
#Merge model averaged results with sagebrush data set so that plot can be generated
glm.ma.sagebrush = cbind(glm.ma.sagebrush, sagebrush)
# plot only the esimates for northern aspect
glm.ma.sagebrush = glm.ma.sagebrush[glm.ma.sagebrush$aspect == "N",]

ggplot(data = glm.ma.sagebrush, aes(x = sagebrush, y = mod.avg.pred)) +
  ggtitle("Naive Occupancy Probability as a function of Sagebrush Density")+
  geom_line()+
  geom_ribbon(aes(ymin = lower.CL, ymax=upper.CL), alpha = .2)+
  xlab("Sagebrush density")+
  ylab("Model averaged naive occupancy probability")+
  ylim(c(0,1))
#dev.off()
##############################################

# predict presence probability for different slope values, while fixing other values
glm.ma.slope <- modavgPred(mods,type="response",newdata=slope,uncond.se="revised")

## manually calculate confidence intervals
logit<-qlogis(glm.ma.slope$mod.avg.pred)
l.se<-glm.ma.slope$uncond.se/(glm.ma.slope$mod.avg.pred*(1-glm.ma.slope$mod.avg.pred))
glm.ma.slope$lower<-plogis(logit-1.96*l.se)
glm.ma.slope$upper<-plogis(logit+1.96*l.se)

#jpeg("PronghornNaiveSlope.jpg",width=800,height=800,res=144)
#Convert model averaged results to a data frame
glm.ma.slope = as.data.frame(glm.ma.slope)
#Merge model averaged results with slope data set so that plot can be generated
glm.ma.slope = cbind(glm.ma.slope, slope)
# plot only the esimates for northern aspect
glm.ma.slope = glm.ma.slope[glm.ma.slope$aspect == "N",]

ggplot(data = glm.ma.slope, aes(x = slope, y = mod.avg.pred)) +
  ggtitle("Naive Occupancy Probability as a function of Percent Slope")+
  geom_line()+
  geom_ribbon(aes(ymin = lower.CL, ymax=upper.CL), alpha = .2)+
  xlab("Percent Slope")+
  ylab("Model averaged naive occupancy probability")+
  ylim(c(0,1))
#dev.off()
##############################################

# predict presence probability for different sagebrush values, while fixing other values
glm.ma.DW <- modavgPred(mods,type="response",newdata=DW,uncond.se="revised")

## manually calculate confidence intervals
logit<-qlogis(glm.ma.DW$mod.avg.pred)
l.se<-glm.ma.DW$uncond.se/(glm.ma.DW$mod.avg.pred*(1-glm.ma.DW$mod.avg.pred))
glm.ma.DW$lower<-plogis(logit-1.96*l.se)
glm.ma.DW$upper<-plogis(logit+1.96*l.se)

#Convert model averaged results to a data frame
glm.ma.DW = as.data.frame(glm.ma.DW)
#Merge model averaged results with distance from water data set so that plot can be generated
glm.ma.DW = cbind(glm.ma.DW, DW)
# plot only the esimates for northern aspect
glm.ma.DW = glm.ma.DW[glm.ma.DW$aspect == "N",]

#jpeg("PronghornNaiveDW.jpg",width=800,height=800,res=144)
ggplot(data = glm.ma.DW, aes(x = DW, y = mod.avg.pred)) +
  ggtitle("Naive Occupancy Probability as a function of Distance from Water")+
  geom_line()+
  geom_ribbon(aes(ymin = lower.CL, ymax=upper.CL), alpha = .2)+
  xlab("Distance from Water")+
  ylab("Model averaged naive occupancy probability")+
  ylim(c(0,1))
#dev.off()
##############################################

# predict presence probability for different sagebrush values, while fixing other values
glm.ma.aspect <- modavgPred(mods,type="response",newdata=aspect,uncond.se="revised")

## manually calculate confidence intervals
logit<-qlogis(glm.ma.aspect$mod.avg.pred)
l.se<-glm.ma.aspect$uncond.se/(glm.ma.aspect$mod.avg.pred*(1-glm.ma.aspect$mod.avg.pred))
glm.ma.aspect$lower<-plogis(logit-1.96*l.se)
glm.ma.aspect$upper<-plogis(logit+1.96*l.se)
#Convert model averaged results to a data frame
glm.ma.aspect = as.data.frame(glm.ma.aspect)
#Merge model averaged results with distance from water data set so that plot can be generated
glm.ma.aspect = cbind(glm.ma.aspect, aspect)

#jpeg("PronghornNaiveAspect.jpg",width=800,height=800,res=144)
ggplot(data = glm.ma.aspect, aes(x = aspect, y = mod.avg.pred)) +
  ggtitle("Naive Occupancy Probability as a function of Aspect")+
  geom_point()+
  geom_errorbar(aes(ymin = lower.CL, ymax=upper.CL), width=.1)+
  xlab("Aspect")+
  ylab("Model averaged naive occupancy probability")+
  ylim(c(0,1))+
  scale_x_discrete(labels=c("North","East","South","West"))
#dev.off()

#######################################################
### use occupancy models to fit 16 models for occupancy component
### with general model for detection component
#######################################################

psi.mod1<-occMod(model=list(psi~1,p~sagebrush+slope+DW+aspect),data=prong,type="so")
psi.mod2<-occMod(model=list(psi~sagebrush,p~sagebrush+slope+DW+aspect),data=prong,type="so")
psi.mod3<-occMod(model=list(psi~slope,p~sagebrush+slope+DW+aspect),data=prong,type="so")
psi.mod4<-occMod(model=list(psi~sagebrush+slope,p~sagebrush+slope+DW+aspect),data=prong,type="so")
psi.mod5<-occMod(model=list(psi~DW,p~sagebrush+slope+DW+aspect),data=prong,type="so")
psi.mod6<-occMod(model=list(psi~sagebrush+DW,p~sagebrush+slope+DW+aspect),data=prong,type="so")
psi.mod7<-occMod(model=list(psi~slope+DW,p~sagebrush+slope+DW+aspect),data=prong,type="so")
psi.mod8<-occMod(model=list(psi~sagebrush+slope+DW,p~sagebrush+slope+DW+aspect),data=prong,type="so")
psi.mod9<-occMod(model=list(psi~aspect,p~sagebrush+slope+DW+aspect),data=prong,type="so")
psi.mod10<-occMod(model=list(psi~sagebrush+aspect,p~sagebrush+slope+DW+aspect),data=prong,type="so")
psi.mod11<-occMod(model=list(psi~slope+aspect,p~sagebrush+slope+DW+aspect),data=prong,type="so")
psi.mod12<-occMod(model=list(psi~sagebrush+slope+aspect,p~sagebrush+slope+DW+aspect),data=prong,type="so")
psi.mod13<-occMod(model=list(psi~DW+aspect,p~sagebrush+slope+DW+aspect),data=prong,type="so")
psi.mod14<-occMod(model=list(psi~sagebrush+DW+aspect,p~sagebrush+slope+DW+aspect),data=prong,type="so")
psi.mod15<-occMod(model=list(psi~slope+DW+aspect,p~sagebrush+slope+DW+aspect),data=prong,type="so")
psi.mod16<-occMod(model=list(psi~sagebrush+slope+DW+aspect,p~sagebrush+slope+DW+aspect),data=prong,type="so")

# compare models by AIC
psi.results<-createAicTable(list(psi.mod1,psi.mod2,psi.mod3,psi.mod4,
                             psi.mod5,psi.mod6,psi.mod7,psi.mod8,
                             psi.mod9,psi.mod10,psi.mod11,psi.mod12,
                             psi.mod13,psi.mod14,psi.mod15,psi.mod16))

# look at AIC table
summary(psi.results)

# regression coefficients for occupancy (psi) and detection (p) from top-ranked model
coef(psi.results$models[[1]], "psi")
coef(psi.results$models[[1]], "p")

#####################################################
# produce some plots of model-averaged estimates to
# illustrate estimated effects
#####################################################

# predict occupancy for different sagebrush values, while fixing other values
# plot only the esimates for northern aspect
psi.ma.sagebrush <- modAvg(psi.results,param="psi",predict=TRUE,newdata=sagebrush[sagebrush$aspect == "N",])

#Merge data sets together for plotting
psi.ma.sagebrush = cbind(psi.ma.sagebrush, sagebrush[sagebrush$aspect=="N",])
#jpeg("PronghornPsiSagebrush.jpg",width=800,height=800,res=144)
ggplot(data = psi.ma.sagebrush, aes(x = sagebrush, y = est)) +
  ggtitle("Model Averaged Occupancy Probability as a function of Sagebrush Density")+
  geom_line()+
  geom_ribbon(aes(ymin = lower_0.95, ymax=upper_0.95), alpha = .2)+
  xlab("Sagebrush density")+
  ylab("Model averaged occupancy probability")+
  ylim(c(0,1))

#dev.off()
##############################################

# predict occupancy for different percent slope values, while fixing other values
# plot only the esimates for northern aspect
psi.ma.slope <- modAvg(psi.results,param="psi",predict=TRUE,newdata=slope[slope$aspect == "N",])

#Merge data sets together for plotting
psi.ma.slope = cbind(psi.ma.slope, slope[slope$aspect=="N",])
#jpeg("PronghornPsislope.jpg",width=800,height=800,res=144)
ggplot(data = psi.ma.slope, aes(x = slope, y = est)) +
  ggtitle("Model Averaged Occupancy Probability as a function of Percent Slope")+
  geom_line()+
  geom_ribbon(aes(ymin = lower_0.95, ymax=upper_0.95), alpha = .2)+
  xlab("Percent Slope")+
  ylab("Model averaged occupancy probability")+
  ylim(c(0,1))
#dev.off()
##############################################
# predict occupancy for different distance to water values, while fixing other values
# plot only the esimates for northern aspect
psi.ma.DW <- modAvg(psi.results,param="psi",predict=TRUE,newdata=DW[DW$aspect == "N",])

#Merge data sets together for plotting
psi.ma.DW = cbind(psi.ma.DW, DW[DW$aspect=="N",])
#jpeg("PronghornPsiDW.jpg",width=800,height=800,res=144)
ggplot(data = psi.ma.DW, aes(x = DW, y = est)) +
  ggtitle("Model Averaged Occupancy Probability as a function of Distance to Water (km)")+
  geom_line()+
  geom_ribbon(aes(ymin = lower_0.95, ymax=upper_0.95), alpha = .2)+
  xlab("Distance to Water (km)")+
  ylab("Model averaged occupancy probability")+
  ylim(c(0,1))
#dev.off()
##############################################
# predict occupancy for different aspects, while fixing other values
psi.ma.aspect <- modAvg(psi.results,param="psi",predict=TRUE,newdata=aspect)

#Create data set for plotting
asp = c("North","East","South","West")
psi.ma.aspect = cbind(psi.ma.aspect, asp)
#jpeg("PronghornPsiaspect.jpg",width=800,height=800,res=144)
ggplot(data = psi.ma.aspect, aes(x = asp, y = est)) +
  ggtitle("Model Averaged Occupancy Probability as a function of Aspect")+
  geom_point()+
  geom_errorbar(aes(ymin = lower_0.95, ymax=upper_0.95), width = .1)+
  xlab("Aspect")+
  ylab("Model averaged occupancy probability")+
  ylim(c(0,1))
#dev.off()

##############################################################
### use occupancy models to fit 16 different models for detection component
### with general model for occupancy component
##############################################################

p.mod1<-occMod(model=list(p~1,psi~sagebrush+slope+DW+aspect),data=prong,type="so")
p.mod2<-occMod(model=list(p~sagebrush,psi~sagebrush+slope+DW+aspect),data=prong,type="so")
p.mod3<-occMod(model=list(p~slope,psi~sagebrush+slope+DW+aspect),data=prong,type="so")
p.mod4<-occMod(model=list(p~sagebrush+slope,psi~sagebrush+slope+DW+aspect),data=prong,type="so")
p.mod5<-occMod(model=list(p~DW,psi~sagebrush+slope+DW+aspect),data=prong,type="so")
p.mod6<-occMod(model=list(p~sagebrush+DW,psi~sagebrush+slope+DW+aspect),data=prong,type="so")
p.mod7<-occMod(model=list(p~slope+DW,psi~sagebrush+slope+DW+aspect),data=prong,type="so")
p.mod8<-occMod(model=list(p~sagebrush+slope+DW,psi~sagebrush+slope+DW+aspect),data=prong,type="so")
p.mod9<-occMod(model=list(p~aspect,psi~sagebrush+slope+DW+aspect),data=prong,type="so")
p.mod10<-occMod(model=list(p~sagebrush+aspect,psi~sagebrush+slope+DW+aspect),data=prong,type="so")
p.mod11<-occMod(model=list(p~slope+aspect,psi~sagebrush+slope+DW+aspect),data=prong,type="so")
p.mod12<-occMod(model=list(p~sagebrush+slope+aspect,psi~sagebrush+slope+DW+aspect),data=prong,type="so")
p.mod13<-occMod(model=list(p~DW+aspect,psi~sagebrush+slope+DW+aspect),data=prong,type="so")
p.mod14<-occMod(model=list(p~sagebrush+DW+aspect,psi~sagebrush+slope+DW+aspect),data=prong,type="so")
p.mod15<-occMod(model=list(p~slope+DW+aspect,psi~sagebrush+slope+DW+aspect),data=prong,type="so")
p.mod16<-occMod(model=list(p~sagebrush+slope+DW+aspect,psi~sagebrush+slope+DW+aspect),data=prong,type="so")

# compare models by AIC
p.results<-createAicTable(list(p.mod1,p.mod2,p.mod3,p.mod4,
                                 p.mod5,p.mod6,p.mod7,p.mod8,
                                 p.mod9,p.mod10,p.mod11,p.mod12,
                                 p.mod13,p.mod14,p.mod15,p.mod16))

# look at AIC table
summary(p.results)

# regression coefficients for occupancy (psi) and detection (p) from top-ranked model
coef(p.results$models[[1]], "psi")
coef(p.results$models[[1]], "p")

#####################################################
# produce some plots of model-averaged estimates to
# illustrate estimated effects
#####################################################

# predict detection for different sagebrush values, while fixing other values
# plot only the esimates for northern aspect
p.ma.sagebrush <- modAvg(p.results,param="p",predict=TRUE,newdata=sagebrush[sagebrush$aspect == "N",])

#Merge data sets together for plotting
p.ma.sagebrush = cbind(p.ma.sagebrush, sagebrush[sagebrush$aspect=="N",])
#jpeg("PronghornPSagebrush.jpg",width=800,height=800,res=144)
ggplot(data = p.ma.sagebrush, aes(x = sagebrush, y = est)) +
  ggtitle("Model Averaged Occupancy Probability as a function of Sagebrush Density")+
  geom_line()+
  geom_ribbon(aes(ymin = lower_0.95, ymax=upper_0.95), alpha = .2)+
  xlab("Sagebrush density")+
  ylab("Model averaged occupancy probability")+
  ylim(c(0,1))

#dev.off()
##############################################
# predict detection for different slope values, while fixing other values
# plot only the esimates for northern aspect
p.ma.slope <- modAvg(p.results,param="p",predict=TRUE,newdata=slope[slope$aspect == "N",])

#Merge data sets together for plotting
p.ma.slope = cbind(p.ma.slope, slope[slope$aspect=="N",])
#jpeg("PronghornPslope.jpg",width=800,height=800,res=144)
ggplot(data = p.ma.slope, aes(x = slope, y = est)) +
  ggtitle("Model Averaged Occupancy Probability as a function of Percent Slope")+
  geom_line()+
  geom_ribbon(aes(ymin = lower_0.95, ymax=upper_0.95), alpha = .2)+
  xlab("Percent Slope")+
  ylab("Model averaged occupancy probability")+
  ylim(c(0,1))
#dev.off()
##############################################
# predict detection for different distance to water values, while fixing other values
# plot only the esimates for northern aspect
p.ma.DW <- modAvg(p.results,param="p",predict=TRUE,newdata=DW[DW$aspect == "N",])

#Merge data sets together for plotting
p.ma.DW = cbind(p.ma.DW, DW[DW$aspect=="N",])
#jpeg("PronghornPDW.jpg",width=800,height=800,res=144)
ggplot(data = p.ma.DW, aes(x = DW, y = est)) +
  ggtitle("Model Averaged Occupancy Probability as a function of Distance to Water (km)")+
  geom_line()+
  geom_ribbon(aes(ymin = lower_0.95, ymax=upper_0.95), alpha = .2)+
  xlab("Distance to Water (km)")+
  ylab("Model averaged occupancy probability")+
  ylim(c(0,1))
#dev.off()
##############################################
# predict detection for different aspects, while fixing other values
p.ma.aspect <- modAvg(p.results,param="p",predict=TRUE,newdata=aspect)

#Create data set for plotting
asp = c("North","East","South","West")
p.ma.aspect = cbind(p.ma.aspect, asp)
#jpeg("PronghornPaspect.jpg",width=800,height=800,res=144)
ggplot(data = p.ma.aspect, aes(x = asp, y = est)) +
  ggtitle("Model Averaged Occupancy Probability as a function of Aspect")+
  geom_point()+
  geom_errorbar(aes(ymin = lower_0.95, ymax=upper_0.95), width = .1)+
  xlab("Aspect")+
  ylab("Model averaged occupancy probability")+
  ylim(c(0,1))
#dev.off()

