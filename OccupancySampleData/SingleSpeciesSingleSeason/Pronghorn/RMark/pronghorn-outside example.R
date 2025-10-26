# load required R libraries

library(RMark)
library(ggplot2)
library(AICcmodavg)

# load in data file
# open this file in excel to examine it's structure if you're not familiar with R  
pronghorn<-read.csv(file.path("..","pronghorn.csv"), header=TRUE, as.is=TRUE, strip.white=TRUE)

# slope covariate in degrees. Convert to proportion slope
pronghorn$slope<-pronghorn$slope/90

# distance to water in m. Convert to km
pronghorn$DW<-pronghorn$DW/1000

#Set aspect as a factor
pronghorn$aspect = as.factor(pronghorn$aspect)

input.history <- data.frame(freq=1,
                            ch=apply(pronghorn[,c("Survey.1","Survey.2")],1,paste, collapse=""), 
                            stringsAsFactors=FALSE)
head(input.history)

input.history = cbind(input.history, pronghorn[,4:7])
head(input.history)

# Change any NA to . in the chapter history
input.history$ch <- gsub("NA",".", input.history$ch, fixed=TRUE)
head(input.history)

#Create the data structure
prong <- process.data(data = input.history, group = "aspect",
                               model = "Occupancy")
summary(prong)

###############################
# prepare some data frames that will be used later for creating graphs 
sagebrush<-expand.grid(sagebrush=seq(0,20,length.out=100),slope=median(pronghorn$slope),DW=median(pronghorn$DW),aspect=c("N","E","S","W"))
slope<-expand.grid(sagebrush=median(pronghorn$sagebrush),slope=seq(0,1,length.out=100),DW=median(pronghorn$DW),aspect=c("N","E","S","W"))
DW<-expand.grid(sagebrush=median(pronghorn$sagebrush),slope=median(pronghorn$slope),DW=seq(0,3,length.out=100),aspect=c("N","E","S","W"))
aspect<-expand.grid(sagebrush=median(pronghorn$sagebrush),slope=median(pronghorn$slope),DW=median(pronghorn$DW),aspect=c("N","E","S","W"))

#################################################################################################
### logistic regression using GLM
### 'naive' analysis that doesn't account for detection probability
#################################################################################################
glm.data = pronghorn[,4:7]

# y =1 if ever detected at plot, =0 if never detected
glm.data$y<-as.numeric(rowSums(pronghorn[,2:3])>0) 

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

#Convert model averaged results to a data frame
glm.ma.slope = as.data.frame(glm.ma.slope)
#Merge model averaged results with slope data set so that plot can be generated
glm.ma.slope = cbind(glm.ma.slope, slope)
# plot only the esimates for northern aspect
glm.ma.slope = glm.ma.slope[glm.ma.slope$aspect == "N",]

#jpeg("PronghornNaiveSlope.jpg",width=800,height=800,res=144)
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

psi.mod1<-RMark::mark(prong,
                      model="Occupancy",
                      model.parameters=list(
                        Psi   =list(formula=~1),
                        p     =list(formula=~sagebrush+slope+DW+aspect) 
                      )
)
psi.mod2<-RMark::mark(prong,
                      model="Occupancy",
                      model.parameters=list(
                        Psi   =list(formula=~sagebrush),
                        p     =list(formula=~sagebrush+slope+DW+aspect) 
                      )
)
psi.mod3<-RMark::mark(prong,
                      model="Occupancy",
                      model.parameters=list(
                        Psi   =list(formula=~slope),
                        p     =list(formula=~sagebrush+slope+DW+aspect) 
                      )
)
psi.mod4<-RMark::mark(prong,
                      model="Occupancy",
                      model.parameters=list(
                        Psi   =list(formula=~sagebrush+slope),
                        p     =list(formula=~sagebrush+slope+DW+aspect) 
                      )
)
psi.mod5<-RMark::mark(prong,
                      model="Occupancy",
                      model.parameters=list(
                        Psi   =list(formula=~DW),
                        p     =list(formula=~sagebrush+slope+DW+aspect) 
                      )
)
psi.mod6<-RMark::mark(prong,
                      model="Occupancy",
                      model.parameters=list(
                        Psi   =list(formula=~sagebrush+DW),
                        p     =list(formula=~sagebrush+slope+DW+aspect) 
                      )
)
psi.mod7<-RMark::mark(prong,
                      model="Occupancy",
                      model.parameters=list(
                        Psi   =list(formula=~slope+DW),
                        p     =list(formula=~sagebrush+slope+DW+aspect) 
                      )
)
psi.mod8<-RMark::mark(prong,
                      model="Occupancy",
                      model.parameters=list(
                        Psi   =list(formula=~sagebrush+slope+DW),
                        p     =list(formula=~sagebrush+slope+DW+aspect) 
                      )
)
psi.mod9<-RMark::mark(prong,
                      model="Occupancy",
                      model.parameters=list(
                        Psi   =list(formula=~aspect),
                        p     =list(formula=~sagebrush+slope+DW+aspect) 
                      )
)
psi.mod10<-RMark::mark(prong,
                       model="Occupancy",
                       model.parameters=list(
                         Psi   =list(formula=~sagebrush+aspect),
                         p     =list(formula=~sagebrush+slope+DW+aspect) 
                       )
)
psi.mod11<-RMark::mark(prong,
                       model="Occupancy",
                       model.parameters=list(
                         Psi   =list(formula=~slope+aspect),
                         p     =list(formula=~sagebrush+slope+DW+aspect) 
                       )
)
psi.mod12<-RMark::mark(prong,
                       model="Occupancy",
                       model.parameters=list(
                         Psi   =list(formula=~sagebrush+slope+aspect),
                         p     =list(formula=~sagebrush+slope+DW+aspect) 
                       )
)
psi.mod13<-RMark::mark(prong,
                       model="Occupancy",
                       model.parameters=list(
                         Psi   =list(formula=~DW+aspect),
                         p     =list(formula=~sagebrush+slope+DW+aspect) 
                       )
)
psi.mod14<-RMark::mark(prong,
                       model="Occupancy",
                       model.parameters=list(
                         Psi   =list(formula=~sagebrush+DW+aspect),
                         p     =list(formula=~sagebrush+slope+DW+aspect) 
                       )
)
psi.mod15<-RMark::mark(prong,
                       model="Occupancy",
                       model.parameters=list(
                         Psi   =list(formula=~slope+DW+aspect),
                         p     =list(formula=~sagebrush+slope+DW+aspect) 
                       )
)
psi.mod16<-RMark::mark(prong,
                       model="Occupancy",
                       model.parameters=list(
                         Psi   =list(formula=~sagebrush+slope+DW+aspect),
                         p     =list(formula=~sagebrush+slope+DW+aspect) 
                       )
)

# compare models by AIC
psi.results<-RMark::collect.models(type = "Occupancy")

# look at AIC table
psi.results

# regression coefficients for occupancy (psi) and detection (p) from top-ranked model
summary(psi.results[[1]])
psi.results[[1]]$results$real

#####################################################
# produce some plots of model-averaged estimates to
# illustrate estimated effects
#####################################################

# predict occupancy for different sagebrush values with North Aspect,
# while fixing other values
Psi.ma <- RMark::model.average(psi.results, param="Psi")
Psi.ma

ddl = make.design.data(prong)
ddl$Psi # see the index numbers

psi.ma.sagebrush <- covariate.predictions(psi.results, indices=10, 
                                          data=sagebrush[sagebrush$aspect=="N",])

#jpeg("PronghornPsiSagebrush.jpg",width=800,height=800,res=144)
ggplot(data = psi.ma.sagebrush$estimates, aes(x = sagebrush, y = estimate)) +
  ggtitle("Model Averaged Occupancy Probability as a function of Sagebrush Density")+
  geom_line()+
  geom_ribbon(aes(ymin = lcl, ymax=ucl), alpha = .2)+
  xlab("Sagebrush density")+
  ylab("Model averaged occupancy probability")+
  ylim(c(0,1))
#dev.off()
##############################################

# predict occupancy for different slope values in the North Aspect,
# while fixing other values
psi.ma.slope <- covariate.predictions(psi.results, indices=10, 
                                          data=slope[slope$aspect=="N",])

#jpeg("PronghornPsislope.jpg",width=800,height=800,res=144)
ggplot(data = psi.ma.slope$estimates, aes(x = slope, y = estimate)) +
  ggtitle("Model Averaged Occupancy Probability as a function of Percent Slope")+
  geom_line()+
  geom_ribbon(aes(ymin = lcl, ymax=ucl), alpha = .2)+
  xlab("Percent Slope")+
  ylab("Model averaged occupancy probability")+
  ylim(c(0,1))
#dev.off()
##############################################

# predict occupancy for different distance to water values in the North Aspect,
# while fixing other values
psi.ma.DW <- covariate.predictions(psi.results, indices=10, 
                                      data=DW[DW$aspect=="N",])

#jpeg("PronghornPsiDW.jpg",width=800,height=800,res=144)
ggplot(data = psi.ma.DW$estimates, aes(x = DW, y = estimate)) +
  ggtitle("Model Averaged Occupancy Probability as a function of Distance to Water (km)")+
  geom_line()+
  geom_ribbon(aes(ymin = lcl, ymax=ucl), alpha = .2)+
  xlab("Distance to Water (km)")+
  ylab("Model averaged occupancy probability")+
  ylim(c(0,1))
#dev.off()
##############################################

# predict occupancy for different aspects, while fixing other values
psi.ma.aspect.E <- covariate.predictions(psi.results, indices=9, 
                                   data=aspect[aspect$aspect == "E",])
psi.ma.aspect.N <- covariate.predictions(psi.results, indices=10, 
                                         data=aspect[aspect$aspect == "N",])
psi.ma.aspect.S <- covariate.predictions(psi.results, indices=11, 
                                         data=aspect[aspect$aspect == "S",])
psi.ma.aspect.W <- covariate.predictions(psi.results, indices=12, 
                                         data=aspect[aspect$aspect == "W",])
psi.ma.aspect = rbind(psi.ma.aspect.E$estimates[,7:11],
                      psi.ma.aspect.N$estimates[,7:11],
                      psi.ma.aspect.S$estimates[,7:11],
                      psi.ma.aspect.W$estimates[,7:11])

#jpeg("PronghornPsiaspect.jpg",width=800,height=800,res=144)
ggplot(data = psi.ma.aspect, aes(x = aspect, y = estimate)) +
  ggtitle("Model Averaged Occupancy Probability as a function of Aspect")+
  geom_point()+
  geom_errorbar(aes(ymin = lcl, ymax=ucl), width = .1)+
  xlab("Aspect")+
  ylab("Model averaged occupancy probability")+
  ylim(c(0,1))+
  scale_x_discrete(labels=c("North","East","South","West"))
#dev.off()

##############################################################
### use occupancy models to fit 16 different models for detection component
### with general model for occupancy component
##############################################################

p.mod1<-RMark::mark(prong,
                      model="Occupancy",
                      model.parameters=list(
                        Psi   =list(formula=~sagebrush+slope+DW+aspect),
                        p     =list(formula=~1) 
                      )
)
p.mod2<-RMark::mark(prong,
                      model="Occupancy",
                      model.parameters=list(
                        Psi   =list(formula=~sagebrush+slope+DW+aspect),
                        p     =list(formula=~sagebrush) 
                      )
)
p.mod3<-RMark::mark(prong,
                      model="Occupancy",
                      model.parameters=list(
                        Psi   =list(formula=~sagebrush+slope+DW+aspect),
                        p     =list(formula=~slope) 
                      )
)
p.mod4<-RMark::mark(prong,
                      model="Occupancy",
                      model.parameters=list(
                        Psi   =list(formula=~sagebrush+slope+DW+aspect),
                        p     =list(formula=~sagebrush+slope) 
                      )
)
p.mod5<-RMark::mark(prong,
                      model="Occupancy",
                      model.parameters=list(
                        Psi   =list(formula=~sagebrush+slope+DW+aspect),
                        p     =list(formula=~DW) 
                      )
)
p.mod6<-RMark::mark(prong,
                      model="Occupancy",
                      model.parameters=list(
                        Psi   =list(formula=~sagebrush+slope+DW+aspect),
                        p     =list(formula=~sagebrush+DW) 
                      )
)
p.mod7<-RMark::mark(prong,
                      model="Occupancy",
                      model.parameters=list(
                        Psi   =list(formula=~sagebrush+slope+DW+aspect),
                        p     =list(formula=~slope+DW) 
                      )
)
p.mod8<-RMark::mark(prong,
                      model="Occupancy",
                      model.parameters=list(
                        Psi   =list(formula=~sagebrush+slope+DW+aspect),
                        p     =list(formula=~sagebrush+slope+DW) 
                      )
)
p.mod9<-RMark::mark(prong,
                      model="Occupancy",
                      model.parameters=list(
                        Psi   =list(formula=~sagebrush+slope+DW+aspect),
                        p     =list(formula=~aspect) 
                      )
)
p.mod10<-RMark::mark(prong,
                       model="Occupancy",
                       model.parameters=list(
                         Psi   =list(formula=~sagebrush+slope+DW+aspect),
                         p     =list(formula=~sagebrush+aspect) 
                       )
)
p.mod11<-RMark::mark(prong,
                       model="Occupancy",
                       model.parameters=list(
                         Psi   =list(formula=~sagebrush+slope+DW+aspect),
                         p     =list(formula=~slope+aspect) 
                       )
)
p.mod12<-RMark::mark(prong,
                       model="Occupancy",
                       model.parameters=list(
                         Psi   =list(formula=~sagebrush+slope+DW+aspect),
                         p     =list(formula=~sagebrush+slope+aspect) 
                       )
)
p.mod13<-RMark::mark(prong,
                       model="Occupancy",
                       model.parameters=list(
                         Psi   =list(formula=~sagebrush+slope+DW+aspect),
                         p     =list(formula=~DW+aspect) 
                       )
)
p.mod14<-RMark::mark(prong,
                       model="Occupancy",
                       model.parameters=list(
                         Psi   =list(formula=~sagebrush+slope+DW+aspect),
                         p     =list(formula=~sagebrush+DW+aspect) 
                       )
)
p.mod15<-RMark::mark(prong,
                       model="Occupancy",
                       model.parameters=list(
                         Psi   =list(formula=~sagebrush+slope+DW+aspect),
                         p     =list(formula=~slope+DW+aspect) 
                       )
)
p.mod16<-RMark::mark(prong,
                       model="Occupancy",
                       model.parameters=list(
                         Psi   =list(formula=~sagebrush+slope+DW+aspect),
                         p     =list(formula=~sagebrush+slope+DW+aspect) 
                       )
)

#Delete old occupancy models before creating new AIC table
rm(psi.mod1,psi.mod2,psi.mod3,psi.mod4,psi.mod5,psi.mod6,psi.mod7,psi.mod8,
   psi.mod9,psi.mod10,psi.mod11,psi.mod12,psi.mod13,psi.mod14,psi.mod15,psi.mod16)

# compare models by AIC
p.results<-RMark::collect.models(type = "Occupancy")

# look at AIC table
p.results

# regression coefficients for occupancy (psi) and detection (p) from top-ranked model
summary(p.results[[1]])
p.results[[1]]$results$real

#####################################################
# produce some plots of model-averaged estimates to
# illustrate estimated effects
#####################################################

# predict detection probability for different sagebrush values with North Aspect,
# while fixing other values
p.ma <- RMark::model.average(p.results, param="p")
p.ma

ddl = make.design.data(prong)
ddl$p # see the index numbers

p.ma.sagebrush <- covariate.predictions(p.results, indices=3, 
                                          data=sagebrush[sagebrush$aspect=="N",])

#jpeg("PronghornPSagebrush.jpg",width=800,height=800,res=144)
ggplot(data = p.ma.sagebrush$estimates, aes(x = sagebrush, y = estimate)) +
  ggtitle("Model Averaged Detection Probability as a function of Sagebrush Density")+
  geom_line()+
  geom_ribbon(aes(ymin = lcl, ymax=ucl), alpha = .2)+
  xlab("Sagebrush density")+
  ylab("Model averaged detection probability")+
  ylim(c(0,1))
#dev.off()
##############################################

# predict detection probability for different slope values in the North Aspect,
# while fixing other values
p.ma.slope <- covariate.predictions(p.results, indices=3, 
                                        data=slope[slope$aspect=="N",])

#jpeg("PronghornPslope.jpg",width=800,height=800,res=144)
ggplot(data = p.ma.slope$estimates, aes(x = slope, y = estimate)) +
  ggtitle("Model Averaged Detection Probability as a function of Percent Slope")+
  geom_line()+
  geom_ribbon(aes(ymin = lcl, ymax=ucl), alpha = .2)+
  xlab("Percent Slope")+
  ylab("Model averaged detection probability")+
  ylim(c(0,1))
#dev.off()
##############################################

# predict detection probability for different distance to water values in the North Aspect,
# while fixing other values
p.ma.DW <- covariate.predictions(p.results, indices=3, 
                                    data=DW[DW$aspect=="N",])

#jpeg("PronghornPDW.jpg",width=800,height=800,res=144)
ggplot(data = p.ma.DW$estimates, aes(x = DW, y = estimate)) +
  ggtitle("Model Averaged Detection Probability as a function of Distance to Water (km)")+
  geom_line()+
  geom_ribbon(aes(ymin = lcl, ymax=ucl), alpha = .2)+
  xlab("Distance to Water (km)")+
  ylab("Model averaged detection probability")+
  ylim(c(0,1))
#dev.off()
##############################################

# predict detection probabilities for different aspects, while fixing other values
p.ma.aspect.E <- covariate.predictions(p.results, indices=1, 
                                         data=aspect[aspect$aspect == "E",])
p.ma.aspect.N <- covariate.predictions(p.results, indices=3, 
                                         data=aspect[aspect$aspect == "N",])
p.ma.aspect.S <- covariate.predictions(p.results, indices=5, 
                                         data=aspect[aspect$aspect == "S",])
p.ma.aspect.W <- covariate.predictions(p.results, indices=7, 
                                         data=aspect[aspect$aspect == "W",])
p.ma.aspect = rbind(p.ma.aspect.E$estimates[,7:11],
                      p.ma.aspect.N$estimates[,7:11],
                      p.ma.aspect.S$estimates[,7:11],
                      p.ma.aspect.W$estimates[,7:11])

#jpeg("PronghornPaspect.jpg",width=800,height=800,res=144)
ggplot(data = p.ma.aspect, aes(x = aspect, y = estimate)) +
  ggtitle("Model Detection Occupancy Probability as a function of Aspect")+
  geom_point()+
  geom_errorbar(aes(ymin = lcl, ymax=ucl), width = .1)+
  xlab("Aspect")+
  ylab("Model averaged detection probability")+
  ylim(c(0,1))+
  scale_x_discrete(labels=c("North","East","South","West"))
#dev.off()

#cleanup
cleanup(ask=FALSE)
