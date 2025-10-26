# Single Species Single Season Occupancy models using JAGS packages

# Blue Gross Beaks.
#Downloaded from https://sites.google.com/site/asrworkshop/home/schedule/r-occupancy-1

#An occupancy study was made on Blue Grosbeaks (Guiraca caerulea) 
# on 41 old fields planted to longleaf pines (Pinus palustris) 
# in southern Georgia, USA. 

# Surveys were 500 m transects across each field 
# and were completed three times during the breeding season in 2001.

# Columns in the file are:
#    field - field number
#    v1, v2, v3 -  detection histories for each site on each of 3 visit during the 2001 breeding season.    
#    field.size - size of the files
#    bqi - Enrollment in bobwihte quail initiative; does occupancy increase if field belongs to this initiative?
#    crop.hist - crop history
#    crop1, crop2 - indicator variables for the crop history
#    count1, count2, count3 - are actual counts of birds detected in each visit


# 2018-08-15 Code contributed by Carl James Schwarz (cschwarz.stat.sfu.cs@gmail.com)


#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#
#    Bayesian model using JAGS with NO COVARIATES, either for psi or detection

library("R2jags")  # used for call to JAGS
library(coda)
library(ggplot2)
library(reshape2)

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
#     psi     - occupancy parameter
#     p       - common detection probability

# 
cat(file="model.txt", "
############################################################

model {
   # set up the state model, i.e. is the site actually occupied or not
   for(i in 1:Nsites){
      z[i] ~  dbern(psi)
   }

   # the observation model.
   for(j in 1:Nsites.visits){
      p.z[j] <- z[Site[j]]*p
      History[j] ~ dbern(p.z[j])
   }

   # priors
   psi ~ dbeta(1,1)
   p   ~ dbeta(1,1)

   # derived variables
   # number of occupied sites
   occ.sites <- sum(z[1:Nsites])
 
   # belief that psi is above some value
   prob.psi.greater.50 <- ifelse( psi > 0.5, 1, 0)
}
") # End of the model



# Next create the data.txt file.
# Initialize the data values using standard R code by either reading
# in from an external file, or plain assignment.
 
input.data <- read.csv(file.path("..","blgr.csv"), 
                       header=TRUE, as.is=TRUE, strip.white=TRUE) 
head(input.data)

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.data)
range(input.data[, c("v1","v2","v3")], na.rm=TRUE)
sum(is.na(input.data[, c("v1","v2","v3")]))

input.history <- input.data[, c("v1","v2","v3")]
head(input.history)

History <- as.vector(unlist(input.history))  # stacks the columns
Site    <- rep(1:nrow(input.history), ncol(input.history))
Visit   <- rep(1:ncol(input.history), each=nrow(input.history))
cbind(Site, Visit, History)[1:10,]

Nsites        <- nrow(input.history)
Nvisits       <- ncol(input.history)
Nsites.visits <- length(History)

# The datalist will be passed to JAGS with the names of the data
# values.
data.list <- list("Nsites","Nvisits","Nsites.visits",
                  "History", "Site", "Visit") # or

# check the list
data.list


# Next create the initial values.
# If you are using more than one chain, you need to create a function
# that returns initial values for each chain.

# We define the initial value of z as 1 if any visit resulted in a detection, other wise 0
init.z <- apply(input.history, 1, max, na.rm=TRUE)

# we will use the naive estimate of occupancy.
init.psi <- sum(init.z)/Nsites

# initial p will be what fraction of 1's exist/occupancy
init.p <- mean(History)/init.psi

# we will start at the same initial starting point for each chain even though this
# is not recommended. 
init.list <- list(
      list(z=init.z, p=init.p, psi=init.psi ),
      list(z=init.z, p=init.p, psi=init.psi ),
      list(z=init.z, p=init.p, psi=init.psi )
      )  # end of list of lists of initial values

# Next create the list of parameters to monitor.
# The deviance is automatically monitored.
# 
monitor.list <- c("z","p", "psi", "occ.sites", "prob.psi.greater.50") # parameters to monitor
 
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
      DIC=TRUE,               # is DIC to be computed?
      working.dir=getwd()    # store results in current working directory
      )




#######################################
# extract some of the usual stuff and use R code directly
# use the standard print method

names(results)
names(results$BUGSoutput)

# get the summary table
results$BUGSoutput$summary
results$BUGSoutput$summary[,c("mean", "sd", "2.5%","97.5%","Rhat", "n.eff")]
results$BUGSoutput$summary[,c("mean", "sd")]


# get just the means
results$BUGSoutput$mean
results$BUGSoutput$mean$psi

# the results$BUGSoutput$sims.array is a 3-d object [iterations, chains, variables]
dim(results$BUGSoutput$sims.array)
results$BUGSoutput$sims.array[1:5,,1:5]
results$BUGSoutput$sims.array[1:5,1,"psi"]


# the results$BUGSoutput$sims.matrix is a 2-d object [iterations, variables] with chains stacked
# on top of each other
dim(results$BUGSoutput$sims.matrix)
results$BUGSoutput$sims.matrix[1:5,1:10]
results$BUGSoutput$sims.matrix[1:5,"psi"]


# make a posterior density plot
plotdata <- data.frame(parm=results$BUGSoutput$sims.matrix[,"psi"], stringsAsFactors=FALSE)
head(plotdata)
postplot.parm <- ggplot2::ggplot( data=plotdata, aes(x=parm, y=..density..))+
  geom_histogram(alpha=0.3)+
  geom_density()+
  ggtitle("Posterior density plot for psi")
postplot.parm
ggsave(plot=postplot.parm, file='psi-posterior-p-dot-psi-dot.png', h=4, w=6, units="in", dpi=300)



# make a trace plot (notice we use the sims.array here)
plotdata <- data.frame(psi=results$BUGSoutput$sims.array[,,"psi"], stringsAsFactors=FALSE)
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
  ggtitle("Trace plot for psi")+
  geom_line(alpha=.2)
traceplot.parm
ggsave(plot=traceplot.parm, file='psi-trace-p-dot-psi-dot.png', h=4, w=6, units="in", dpi=300)


# autocorrelation plot
# First compute the autocorrelation plot
acf.parm <-acf( results$BUGSoutput$sims.matrix[,"psi"], plot=FALSE)
acf.parm
acfplot.parm <- ggplot(data=with(acf.parm, data.frame(lag, acf)), aes(x = lag, y = acf)) +
  ggtitle("Autocorrelation plot for psi")+
  geom_hline(aes(yintercept = 0)) +
  geom_segment(aes(xend = lag, yend = 0))
acfplot.parm
ggsave(plot=acfplot.parm, file="psi-acf-p-dot-psi-dot.png",h=4, w=6, units="in", dpi=300)




#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#
#    Bayesian model using JAGS with COVARIATES for detection
#    such a time effects.
#    You need to build the covariate matrix and send to JAGS

library("R2jags")  # used for call to JAGS
library(coda)
library(ggplot2)
library(reshape2)

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
#     psi     - occupancy parameter
#     p       - detection probability
#
#     Ncovar.p - number of covariates for logit(p)
#     covar.p - matrix of covariates (Nsites.visits x N.covar.p)

# 
cat(file="model.txt", "
############################################################

model {
   # set up the state model, i.e. is the site actually occupied or not
   for(i in 1:Nsites){
      z[i] ~  dbern(psi)
   }

   # the observation model.
   for(j in 1:Nsites.visits){
      logit(p.detect[j]) <- inprod(beta.p[1:Ncovar.p], covar.p[j, 1:Ncovar.p])
      p.z[j] <- z[Site[j]]*p.detect[j]
      History[j] ~ dbern(p.z[j])
   }

   # priors
   psi ~ dbeta(1,1)
   for(i in 1:Ncovar.p){
      beta.p[i]   ~ dnorm(0, .0001)
   }

   # derived variables
   # number of occupied sites
   occ.sites <- sum(z[1:Nsites])
 
   # belief that psi is above some value
   prob.psi.greater.50 <- ifelse( psi > 0.5, 1, 0)
}
") # End of the model



# Next create the data.txt file.
# Initialize the data values using standard R code by either reading
# in from an external file, or plain assignment.
 
input.data <- read.csv(file.path("..","blgr.csv"), 
                       header=TRUE, as.is=TRUE, strip.white=TRUE) 
head(input.data)

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.data)
range(input.data[, c("v1","v2","v3")], na.rm=TRUE)
sum(is.na(input.data[, c("v1","v2","v3")]))

input.history <- input.data[, c("v1","v2","v3")]
head(input.history)

History <- as.vector(unlist(input.history))  # stacks the columns
Site    <- rep(1:nrow(input.history), ncol(input.history))
Visit   <- rep(1:ncol(input.history), each=nrow(input.history))

covar.p <- model.matrix( ~as.factor(Visit), data=as.data.frame(Visit))
covar.p[1:10,]

# create a covariate matrix for time effects.
cbind(Site, Visit, History, covar.p)[1:10,]

Nsites        <- nrow(input.history)
Nvisits       <- ncol(input.history)
Nsites.visits <- length(History)
Ncovar.p      <- ncol(covar.p)

# The datalist will be passed to JAGS with the names of the data
# values.
data.list <- list("Nsites","Nvisits","Nsites.visits",
                  "History", "Site", "Visit",
                  "Ncovar.p","covar.p") # or

# check the list
data.list


# Next create the initial values.
# If you are using more than one chain, you need to create a function
# that returns initial values for each chain.

# We define the initial value of z as 1 if any visit resulted in a detection, other wise 0
init.z <- apply(input.history, 1, max, na.rm=TRUE)

# we will use the naive estimate of occupancy.
init.psi <- sum(init.z)/Nsites

# initial p will be what fraction of 1's exist/occupancy
init.p <- mean(History)/init.psi

# we will start at the same initial starting point for each chain even though this
# is not recommended. 
init.list <- list(
      list(z=init.z, beta.p=rep(0,ncol(input.history)), psi=init.psi ),
      list(z=init.z, beta.p=rep(0,ncol(input.history)), psi=init.psi ),
      list(z=init.z, beta.p=rep(0,ncol(input.history)), psi=init.psi )
      )  # end of list of lists of initial values

# Next create the list of parameters to monitor.
# The deviance is automatically monitored.
# 
monitor.list <- c("z",
                  "p.detect", "beta.p",
                  "psi", "occ.sites", "prob.psi.greater.50") # parameters to monitor
 
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
      DIC=TRUE,               # is DIC to be computed?
      working.dir=getwd()    # store results in current working directory
      )




#######################################
# extract some of the usual stuff and use R code directly
# use the standard print method

names(results)
names(results$BUGSoutput)

# get the summary table
results$BUGSoutput$summary
results$BUGSoutput$summary[,c("mean", "sd", "2.5%","97.5%","Rhat", "n.eff")]
results$BUGSoutput$summary[,c("mean", "sd")]


# get just the means
results$BUGSoutput$mean
results$BUGSoutput$mean$psi

# the results$BUGSoutput$sims.array is a 3-d object [iterations, chains, variables]
dim(results$BUGSoutput$sims.array)
results$BUGSoutput$sims.array[1:5,,1:5]
results$BUGSoutput$sims.array[1:5,1,"psi"]


# the results$BUGSoutput$sims.matrix is a 2-d object [iterations, variables] with chains stacked
# on top of each other
dim(results$BUGSoutput$sims.matrix)
results$BUGSoutput$sims.matrix[1:5,1:10]
results$BUGSoutput$sims.matrix[1:5,"psi"]


# make a posterior density plot
plotdata <- data.frame(parm=results$BUGSoutput$sims.matrix[,"psi"], stringsAsFactors=FALSE)
head(plotdata)
postplot.parm <- ggplot2::ggplot( data=plotdata, aes(x=parm, y=..density..))+
  geom_histogram(alpha=0.3)+
  geom_density()+
  ggtitle("Posterior density plot for psi")
postplot.parm
ggsave(plot=postplot.parm, file='psi-posterior-p-t-psi-dot.png', h=4, w=6, units="in", dpi=300)



# make a trace plot (notice we use the sims.array here)
plotdata <- data.frame(psi=results$BUGSoutput$sims.array[,,"psi"], stringsAsFactors=FALSE)
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
  ggtitle("Trace plot for psi")+
  geom_line(alpha=.2)
traceplot.parm
ggsave(plot=traceplot.parm, file='psi-trace-p-t-psi-dot.png', h=4, w=6, units="in", dpi=300)


# autocorrelation plot
# First compute the autocorrelation plot
acf.parm <-acf( results$BUGSoutput$sims.matrix[,"psi"], plot=FALSE)
acf.parm
acfplot.parm <- ggplot(data=with(acf.parm, data.frame(lag, acf)), aes(x = lag, y = acf)) +
  ggtitle("Autocorrelation plot for psi")+
  geom_hline(aes(yintercept = 0)) +
  geom_segment(aes(xend = lag, yend = 0))
acfplot.parm
ggsave(plot=acfplot.parm, file="psi-acf-p-t-psi-dot.png",h=4, w=6, units="in", dpi=300)

