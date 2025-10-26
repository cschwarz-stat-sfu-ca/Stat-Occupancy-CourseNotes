# Single Species Single Season Occupancy models 

# Weta Data

#  using JAGS

# 2018-08-15 Code contributed by Carl James Schwarz (cschwarz.stat.sfu.cs@gmail.com)

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

library("R2jags")  # used for call to JAGS
library(car)
library(coda)
library(ggplot2)
library(readxl)
library(reshape2)

options(width=400)  # make html output wider

# get the data read in
# Data for detections should be a data frame with rows corresponding to sites
# and columns to visits.
# The usual 1=detected; 0=not detected; NA=not visited is used.

input.history <- readxl::read_excel(file.path("..","weta.xls"),
                                    sheet="detection_histories",
                                    na="-",
                                    col_names=FALSE)  # notice no column names in row 1 of data file. 

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.history)
ncol(input.history)
range(input.history, na.rm=TRUE)
sum(is.na(input.history))

head(input.history)

# Get the site level covariates
site_covar <- readxl::read_excel(file.path("..","weta.xls"),
                                 sheet="site_covar",
                                 na="-",
                                 col_names=TRUE)  # notice col_names in row 1 of table. 


# Create an alternate site level covariate that is a categorical variable rather 
# than indicator variables
site_covar$BrowCat <- paste(c("","B")[1+unlist(site_covar[,1])], c("","N")[1+unlist(site_covar[,2])], sep="")
xtabs(~BrowCat, data=site_covar,exclude=NULL, na.action=na.pass)
colSums(site_covar[,1:2])
site_covar$Site <- 1:nrow(site_covar)

head(site_covar)


# Get the individual covariates. 
obs1 <- readxl::read_excel(file.path("..","weta.xls"),
                           sheet="Obs1",
                           na="-",
                           col_names=FALSE) 
obs2 <- readxl::read_excel(file.path("..","weta.xls"),
                           sheet="Obs2",
                           na="-",
                           col_names=FALSE) 
obs3 <- readxl::read_excel(file.path("..","weta.xls"),
                           sheet="Obs3",
                           na="-",
                           col_names=FALSE) 

Obs <- obs1*1 + obs2*2 + obs3*3
head(Obs)

# Observational covariate needs to be "stacked" so that sites1...siteS for survey occastion 1
# are then followed by covariate at survey occastion 2 for sites1...siteS, etc

survey.cov <- data.frame(site=rep(1:nrow(input.history) , ncol(input.history)),
                         visit=as.character(rep(1:ncol(input.history), each=nrow(input.history))),  # notice we make a character 
                         obs1 =as.vector(unlist(obs1)),
                         obs2 =as.vector(unlist(obs2)),
                         obs3 =as.vector(unlist(obs3)),
                         Obs  =as.character(as.vector(unlist(Obs))),    # notice we make a character string
                         stringsAsFactors=FALSE)
head(survey.cov)
str(survey.cov)

# check that missing values in history and observer covariates align
select <- is.na(as.vector(unlist(input.history)))
survey.cov[select,]
sum(is.na(survey.cov[!select,]))

# The missing values in the survey covariates must be filled with dummy
# values to avoid problems in fitting the models that depend on them
survey.cov[ is.na(survey.cov)] <- -1
survey.cov[select,]
sum(is.na(survey.cov[!select,]==-1))



# The BUGS model is specified as a text file.

# The model file.
# The cat() command is used to save the model to the working directory.
# Notice that you CANNOT have any " (double quotes) in the bugs code
# between the start and end of the cat("...",) command.

# Inputs to the model are 
#     Nsites  - number of sites
#     Nvisits - (max) number of visits over all sites.
#     Nsites.visits - number of sites x number of visits 
#          if there is missing data (no visits), simply drop the corresponding row
#     History - vector of 1 or 0 corresponding to Site-Visit pair
#     Site    - vector indicating which site the row corresponds to
#     Visit   - vector indicating which visit the row corresponds to
# 
#     dmatrix.psi - design matrix for psi
#     Nbeta.psi     - number of columns in design matrix for psi

#     dmatrix.p   - design matrix for p
#     Nbeta.p       - number of columns of design matrix for p

# 
cat(file="model.txt", "
############################################################

model {
   # estimate psi for each site from the design matrix
   for(i in 1:Nsites){
      logit(psi[i]) = inprod( dmatrix.psi[i, 1:Nbeta.psi], beta.psi[1:Nbeta.psi])
   }

   # estimate p for each observation
   for(i in 1:Nsites.visits){
      logit(p[i]) = inprod( dmatrix.p[i, 1:Nbeta.p], beta.p[1:Nbeta.p])
    }


   # set up the state model, i.e. is the site actually occupied or not
   for(i in 1:Nsites){
      z[i] ~  dbern(psi[i])
   }

   # the observation model.
   for(j in 1:Nsites.visits){
      p.z[j] <- z[Site[j]]*p[j]
      History[j] ~ dbern(p.z[j])
   }

   # priors on the betas
   beta.psi[1] ~ dnorm(0, .25)  # intercept we want basically flat on regular scale
   for(i in 2:Nbeta.psi){
      beta.psi[i] ~ dnorm(0, .0001)
   }

   beta.p[1] ~ dnorm(0, .25)
   for(i in 2:Nbeta.p){
      beta.p[i] ~ dnorm(0, .0001)
   }

   # derived variables
   # number of occupied sites
   occ.sites <- sum(z[1:Nsites])
 
   # belief that psi is above some value
   prob.psi.greater.50 <- ifelse( psi[1] > 0.5, 1, 0)
}
") # End of the model


# get the data in the right format. We want all of the sites for visit 1, then all of the sites for visit 2 etc.

Survey <- data.frame(
     History = as.vector(unlist(input.history)),  # stacks the columns
     Site    = rep(1:nrow(input.history), ncol(input.history)),
     Visit   = as.character(rep(1:ncol(input.history), each=nrow(input.history))), 
     stringsAsFactors=FALSE)
Survey[1:10,]

# add in covaraites
Survey$Obs <- survey.cov$Obs  # both are in correct order
Survey <- merge(Survey, site_covar[,c("Site","BrowCat")])

# re-sort
Survey <- Survey[ order(Survey$Visit, Survey$Site),]
head(Survey)
str(Survey) # be sure that all categorical variables are character or factors

# Remove any rows with missing history value (i.e missing)
sum(is.na(Survey$History))
dim(Survey)
Survey <- Survey[!is.na(Survey$History),]
dim(Survey)

Nsites        <- nrow(input.history)
Nvisits       <- ncol(input.history)
Nsites.visits <- nrow(Survey)



# Get the design matrix for model psi(BrowCat)p(visit)
dmatrix.psi <- model.matrix(~BrowCat,  data=site_covar)  
Nbeta.psi <- ncol(dmatrix.psi)

dmatrix.p  <- model.matrix(~Visit,  data=Survey)
Nbeta.p    <- ncol(dmatrix.p)
model.name <- "psi(B) p(v)"


# Get the design matrix for model psi(BrowCat) p(Obs + Visit)
dmatrix.psi <- model.matrix(~BrowCat,  data=site_covar)  # constant psi
Nbeta.psi <- ncol(dmatrix.psi)

dmatrix.p  <- model.matrix(~Visit+Obs,  data=Survey)
Nbeta.p    <- ncol(dmatrix.p)
model.name <- "pt-psidot"







# The datalist will be passed to JAGS with the names of the data
# values.
data.list <- list(Nsites=Nsites,
                  Nvisit=Nvisits,
                  Nsites.visits=Nsites.visits,
                  History = Survey$History,
                  Site    = Survey$Site,
                  Visit   = Survey$Visit,
                  dmatrix.psi=dmatrix.psi, Nbeta.psi=Nbeta.psi,
                  dmatrix.p  =dmatrix.p,   Nbeta.p  =Nbeta.p)
                  
  
# check the list
data.list


# Next create the initial values.
# If you are using more than one chain, you need to create a function
# that returns initial values for each chain.

# We define the initial value of z as 1 if any visit resulted in a detection, other wise 0
init.z <- apply(input.history, 1, max, na.rm=TRUE)

# we will start at the same initial starting point for each chain even though this
# is not recommended. 
init.list <- list(
      list(z=init.z, beta.psi=rep(0,Nbeta.psi), beta.p=rep(0,Nbeta.p) ),
      list(z=init.z, beta.psi=rep(0,Nbeta.psi), beta.p=rep(0,Nbeta.p) ),
      list(z=init.z, beta.psi=rep(0,Nbeta.psi), beta.p=rep(0,Nbeta.p) )
      
)  # end of list of lists of initial values

# Next create the list of parameters to monitor.
# The deviance is automatically monitored.
# 
monitor.list <- c("z","occ.sites", "prob.psi.greater.50",
                  "psi", "p",
                  "beta.psi", "beta.p") # parameters to monitor
 
# Finally, the actual call to JAGS
set.seed(234234)  # intitalize seed for MCMC 

results <- R2jags::jags( 
      data      =data.list,   # list of data variables
      inits     =init.list,   # list/function for initial values
      parameters=monitor.list,# list of parameters to monitor
      model.file="model.txt",  # file with bugs model
      n.chains=3,
      n.iter  =5000,          # total iterations INCLUDING burn in
      n.burnin=2000,          # number of burning iterations
      n.thin=2,               # how much to thin
      DIC=TRUE               # is DIC to be computed?
      )




#######################################
# extract some of the usual stuff and use R code directly
# use the standard print method

names(results)
names(results$BUGSoutput)

# get the summary table
results$BUGSoutput$summary
#results$BUGSoutput$summary[,c("mean", "sd", "2.5%","97.5%","Rhat", "n.eff")]
#results$BUGSoutput$summary[,c("mean", "sd")]


# get just the means
results$BUGSoutput$mean
results$BUGSoutput$mean$psi

# the results$BUGSoutput$sims.array is a 3-d object [iterations, chains, variables]
dim(results$BUGSoutput$sims.array)
results$BUGSoutput$sims.array[1:5,1,1:10]
results$BUGSoutput$sims.array[1:5,1,"psi[1]", drop=FALSE]


# the results$BUGSoutput$sims.matrix is a 2-d object [iterations, variables] with chains stacked
# on top of each other
dim(results$BUGSoutput$sims.matrix)
results$BUGSoutput$sims.matrix[1:5,1:10]
results$BUGSoutput$sims.matrix[1:5,"psi[1]", drop=FALSE]


# make a posterior density plot
plotdata <- data.frame(parm=results$BUGSoutput$sims.matrix[,c("psi[1]","psi[4]")], stringsAsFactors=FALSE) # browse and unbrowsed
head(plotdata)
plotdata2 <- reshape2::melt(plotdata, variable.name="Site", value.name="prob")
head(plotdata2)
postplot.parm <- ggplot2::ggplot( data=plotdata2, aes(x=prob, y=..density..))+
  geom_histogram(alpha=0.3)+
  geom_density()+
  ggtitle("Posterior density plot for psi[1] and psi[[4]")+
  facet_wrap(~Site, ncol=1)
postplot.parm
ggsave(plot=postplot.parm, 
       file=paste('psi-posterior-',model.name,'.png',sep=""), h=4, w=6, units="in", dpi=300)


# Get the odds ratio for Browser vs Not Browse
plotdata <- data.frame(parm=results$BUGSoutput$sims.matrix[,c("beta.psi[2]")], stringsAsFactors=FALSE) # browse and unbrowsed
head(plotdata)
plotdata$odds.ratio <- exp(plotdata$parm)
range(plotdata$odds.ratio)  # some very large odds ratios
summary(plotdata$odds.ratio)
quantile(plotdata$odds.ratio, prob=.98)
plotdata$odds.ratio[ plotdata$odds.ratio > 5] <- NA
head(plotdata)
oddsplot.parm <- ggplot2::ggplot( data=plotdata, aes(x=odds.ratio, y=..density..))+
  geom_histogram(alpha=0.3)+
  geom_density()+
  ggtitle("Odds ratio of occupancy(not browsed):occupancy(browsed)]")
oddsplot.parm
ggsave(plot=postplot.parm, 
       file=paste('odds-psi-posterior-',model.name,'.png',sep=""), h=4, w=6, units="in", dpi=300)



# make a trace plot (notice we use the sims.array here)
plotdata <- data.frame(psi=results$BUGSoutput$sims.array[,,"psi[1]"], stringsAsFactors=FALSE)
plotdata$iteration <- 1:nrow(plotdata)
head(plotdata)

# convert from wide to long format
plotdata2 <- reshape2:::melt.data.frame(data=plotdata, 
                            id.vars="iteration",
                            measure.vars=paste("psi",1:results$BUGSoutput$n.chains,sep="."),
                            variable.name="chain",
                            value.name='psi')
head(plotdata2)
traceplot.parm <- ggplot2::ggplot(data=plotdata2, aes(x=iteration, y=psi, color=chain))+
  ggtitle("Trace plot for psi[1]")+
  geom_line(alpha=.2)
traceplot.parm
ggsave(plot=traceplot.parm, 
       file=paste('psi-trace-',model.name,'.png',sep=""), h=4, w=6, units="in", dpi=300)


# autocorrelation plot
# First compute the autocorrelation plot
acf.parm <-acf( results$BUGSoutput$sims.matrix[,"psi[1]"], plot=FALSE)
acf.parm
acfplot.parm <- ggplot(data=with(acf.parm, data.frame(lag, acf)), aes(x = lag, y = acf)) +
  ggtitle("Autocorrelation plot for psi[1]")+
  geom_hline(aes(yintercept = 0)) +
  geom_segment(aes(xend = lag, yend = 0))
acfplot.parm
ggsave(plot=acfplot.parm, 
       file=paste("psi-acf-", model.name,".png",sep=""),h=4, w=6, units="in", dpi=300)

