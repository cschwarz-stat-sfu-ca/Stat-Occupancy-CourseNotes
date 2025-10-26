# multi season multi-species model

# Fidino M, Simonis JL, Magle SB (2018) 
# A multi-species dynamic occupancy model to estimate local colonization-extinction 
# rates and patterns of co-occurrence between two or more interacting species. 
# Methods in Ecology and Evolution, online in advance of print. 
# https://doi.org/10.1111/2041-210x.13117
# 
# Several years of data on occurrences of 3 species. See their paper.
# 
# Data extracted from Dryad Digital Repository. 
# https://doi.org/10.5061/dryad.k5mp137

# 2018-12-16 Code contributed by Carl James Schwarz (cschwarz.stat.sfu.cs@gmail.com)


#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------

library(car)
library(ggplot2)
library(readxl)
library(reshape2)
library(RMark)



# Robust design, Multi species, Reproduction parameterization
setup.parameters("RDMSOccRepro", check=TRUE) # returns a vector of parameter names (case sensitive)



# phi0  = initial occupancy probabilities in each state parameterized as  phi0[1]*(1-phi0[2]), phi0[1](phi0[2]])

#    1 - Phi0[1]
#   Phi0[1](1 - Phi0[2])
#   Phi0[1]*Phi0[2]

# psi   = p(occupied in year t+1 | not occupied, occupied but not breeding, occupied and breeding in year t)
# r   = p(reproduction in year t+1 | not occupied, occupied but not breeding, occupied and breeding in year t+1)
# p   = p(detection in year t | state 1, state 2 in yar t
# delta=p(detect in state 2 | in state 2)

# The parameters Psi[0], Psi[1], Psi[2], R[0,2], R[1,2], and R[2,2] are used to construct 
# the Phi_t matrix for transitions each interval (midway down second column page 826 of MacKenzie et al. 2009):
#  
#  1 - Psi[0](t)          Psi[0](t)*{1 - R[0,2](t)}          Psi[0](t)*R[0,2](t)
#  1 - Psi[1](t)          Psi[1](t)*{1 - R[1.2](t)}          Psi[1](t)*R[1,2](t)
#  1 - Psi[2](t)          Psi[2](t)*{1 - R[2,2](t)}          Psi[2](t)*R[2,2](t)

# These 3 PIMs for each primary interval are used to construct the detection matrix 
# (top of first column page 826 MacKenzie et al. 2009):
  
#                                                    Observed State                                        
#                                    0                    1                           2 
#------------------------------------------------------------------------------------------------
#True State           0|             1                    0                           0         
#                     1|             1 - p[1]            p[1]                         0         
#                     2|             1 - p[2]            p[2]*(1 - Delta[2,2])    p[2]*Delta[2,2]


input.history <- read.csv(file.path("..","Fidino_MEE","data","fidino_sp_data.csv"),
                          as.is=TRUE, strip.white=TRUE, header=TRUE)
head(input.history)

hawk.history <- data.frame(ch        =apply(input.history[,-(1:2)],1,paste,sep="", collapse=""),
                        freq      =1,
                        tsh       =input.history$tsh,
                        Territory = input.history$TerritoryN,
                        stringsAsFactors=FALSE)
head(hawk.history)

# fix up the NA in the history
hawk.history$ch <- gsub("NA",".", hawk.history$ch)
head(hawk.history)

# make a plot of the imputed occupancy
plotdata <- reshape2::melt(input.history,
                           id.vars=c("TerritoryN","tsh"),
                           value.name="ObsOcc",
                           variable.name="Visit")
plotdata$Visit <- as.numeric(substring(plotdata$Visit,2))
head(plotdata)

ggplot(data=plotdata, aes(x=Visit, y=TerritoryN, color=as.factor(ObsOcc), shape=as.factor(ObsOcc)))+
    ggtitle("Observed occupancy state", subtitle="Not corrected for false negatives")+
    geom_point()+
    scale_shape_manual(name="Observed\nOccupancy", values=c(1, 16, 19))+
    scale_color_manual(name="Observed\nOccupancy", values=c("black","green","red"))+
    theme_bw()
 

hawk.data <- process.data(data=hawk.history,
                          model="RDMSOccRepro",
                          time.intervals=c( rep( c(rep(0,8),1),5),rep(0,8))) # 6 years x 9 visit occasion/year

summary(hawk.data)

# add visit level covariates
# In this case none
hawk.ddl <- make.design.data(hawk.data)
hawk.ddl

# Robust design, Multi State, Reproduction parameterization
setup.parameters("RDMSOccRepro", check=TRUE) # returns a vector of parameter names (case sensitive)

model.list.csv <- textConnection("
Phi0,            Psi,           R,            p,              Delta  
~stratum,       ~stratum,       ~stratum,     ~stratum,         ~1
~stratum*tsh,   ~stratum,       ~stratum,     ~stratum,         ~1
~stratum*tsh,   ~stratum*tsh,   ~stratum*tsh, ~stratum,         ~1
~stratum*tsh,   ~stratum,       ~1,           ~stratum,         ~1
~stratum,       ~stratum*tsh,   ~stratum,     ~stratum,         ~1
~stratum*tsh,   ~stratum*tsh,   ~stratum,     ~stratum,         ~1") 

model.list <- read.csv(model.list.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)
model.list$model.number <- 1:nrow(model.list)
model.list

# fit the models
myobs <- ls()
myobs <- myobs[ grepl("m...",myobs,fixed=TRUE)]
cat("Removing ", myobs, "\n")
rm(list=myobs)

model.fits <- plyr::dlply(model.list, "model.number", function(x,hawk.data,input.ddl){
    cat("\n\n***** Starting ", unlist(x), "\n")
    #browser()
  
      #fit <- myoccMod(model=list(as.formula(paste("psi",x$psi)),
    fit <- RMark::mark(hawk.data, ddl = input.ddl,
                       model="RDMSOccRepro",
                       model.parameters=list(
                           Phi0  =list(formula=as.formula(eval(x$Phi0))),
                           Psi   =list(formula=as.formula(eval(x$Psi))),
                           R     =list(formula=as.formula(eval(x$R))),
                           p     =list(formula=as.formula(eval(x$p))),
                           Delta =list(formula=as.formula(eval(x$Delta))))
                           ,brief=TRUE,output=FALSE, delete=TRUE
                          #,invisible=TRUE,output=TRUE  # set for debugging
            )
    mnumber <- paste("m...",formatC(x$model.number, width = 3, format = "d", flag = "0"),sep="")
    assign( mnumber, fit, envir=.GlobalEnv)
      #browser()
    fit

},hawk.data=hawk.data, input.ddl = hawk.ddl)


model.number <-2

summary(model.fits[[model.number]])
model.fits[[model.number]]$results$real
model.fits[[model.number]]$results$beta


model.fits[[model.number]]$results$derived

get.real(model.fits[[model.number]], "Phi0", se=TRUE)
get.real(model.fits[[model.number]], "Psi", se=TRUE)
get.real(model.fits[[model.number]], "R",    se=TRUE)
temp <- get.real(model.fits[[model.number]], "p",    se=TRUE)
temp[ temp$time==1 & temp$session==1,]

# collect models and make AICc table

model.set <- RMark::collect.models( type="RDMSOccRepro")
model.set

# get model averaged values for the initial occupancy probabilities at different values of the covariates
get.real(model.fits[[model.number]], "Phi0", se=TRUE)


Phi0.ma <- RMark::model.average(model.set, param="Phi0")
Phi0.ma

range(hawk.history$tsh, na.rm=TRUE)
Phi0.pred <- covariate.predictions(model.set, data=hawk.history, indices=1:4)
Phi0.pred$estimates[1:10,]

# Extract the estimated psi values and combine back with standardized values
# This is tricky because the region creates different group
plotdata1 <- Phi0.pred$estimates[ Phi0.pred$estimates$par.index==1, ]
plotdata1$Source <- "Initial occupancy (breeding or not breeding)"

plotdata2 <-  Phi0.pred$estimates[ Phi0.pred$estimates$par.index==2, ]
plotdata2$Source <- "Conditional p(breeding | occupied)"

plotdata <- rbind(plotdata1, plotdata2)
plotdata.ci <- plyr::ddply(plotdata, c("Source"), function(x){
  # select 5 points at random for ci
  x <- x[ sample(nrow(x),5),]
  x
})

head(plotdata)
i.occ.plot <- ggplot(data=plotdata, aes(x=tsh, y=estimate))+
  ggtitle(paste('Model averaged initial occupancy probability from multi-state occupancy model',
                 '\nVariation in trend represents effects from other models',sep=""))+
  geom_point()+
  ylim(0,1)+
  ylab("Initial occupancy probability\nSelected 95% ci")+
  xlab("Area of suitable habitat ")+
  theme(legend.position=c(1,0), legend.justification=c(1,0))+
  facet_wrap(~Source, ncol=1)+
  geom_errorbar(data=plotdata.ci, aes(ymin=lcl, ymax=ucl), width=0.1)
i.occ.plot


# Model averaged estimates of occupancy over time given previous state
Psi.ma <- RMark::model.average(model.set, param="Psi")
Psi.ma


R.ma <- RMark::model.average(model.set, param="R")
R.ma


# Get the transition matrix
model.number <- 2
names(model.fits[[model.number]]$results$derived)
model.fits[[model.number]]$results$real

model.fits[[model.number]]$results$derived$'Phi_t transition matrix'

# Model average the derived parameter
source(file.path("..","..","..","AdditionalFunctions","RMark.additional.functions.R"))

trans.ma <- RMark.model.average.derived(model.set, "Phi_t transition matrix")
head(trans.ma)

temp <- matrix(trans.ma[1:9,"estimate"], 3,3, byrow=TRUE)
temp


# Estimate the average level 1 occupancy over time
# This is a derived parameter and so we need to model average as before.

psi1.ma <- RMark.model.average.derived(model.set, parameter="psi")
psi1.ma$Year <- 1:6
psi1.ma$Source="Model"
psi1.ma$Level="Level 1 occupancy"
psi1.ma

psi2.ma <- RMark.model.average.derived(model.set, parameter="psi_2")
psi2.ma$Year <- 1:6
psi2.ma$Source="Model"
psi2.ma$Level="Level 2 occupancy (Breeding)"
psi2.ma


plotdata<- plyr::rbind.fill(psi1.ma,  psi2.ma )
plotdata

plot.occ1 <- ggplot(data=plotdata, aes(x=Year, y=estimate, linetype=Source))+
  ggtitle("Model averaged level 1 and 2 occupancy from multi-state model", subtitle="At a global average tsh")+
  geom_line()+
  ylim(0,1)+ylab("Average level 1 occupancy over all territories \n(95% ci)")+
  theme(legend.position=c(0,.5), legend.justification=c(0,1), legend.box="horizontal")+
  scale_linetype_discrete(name="Estimate")+
  geom_errorbar(aes(ymin=pmax(0,lcl), ymax=pmin(1,ucl),width=.1))+
  facet_wrap(~Level, ncol=1)
plot.occ1




remove.mark(collect.models(), 1:nrow(collect.models()$model.table))

cleanup(ask=FALSE)
