#    Single Species Multi Season   

# Northern Spotted Owl (Strix occidentalis caurina) in California.
# s=55$ sites visited up to K=$ times per season between 1997 and 2001 ($Y=5$).
# Detection probabilities relatively constant within years, but likely different among years.

# Using JAGS

# 2018-08-17 Code contributed by Carl James Schwarz (cschwarz.stat.sfu.ca@gmail.com)
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------


#    Bayesian model using JAGS with COVARIATES for all effects
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
#     Nyears  - number of years (seasons)
#     Nvisits - (max) number of visits over all sites.

#     Nsites.visits - number of sites x number of visits (overall all years)
#          if there is missing data (no visits), simply drop the corresponding row
#     History - vector of 1 or 0 corresponding to Site-Visit pair
#     Site    - vector indicating which site the history value  corresponds to
#     Year    - vector indicating which year the history value  corresponds to
#     Visit   - vector indicating which visit the history value corresponds to
# 
#     psi     - occupancy parameter
#     p       - detection probability
#
#     *** The initial psi covariate matrix has 1 initial column with Site #
#     Ncovar.psi - number of covariates for logit(psi)
#     covar.psi  - matrix of covariates (Nsites x N.covar.psi)

#     *** ALL of the following covariate matrices have two initial columns with Site and Year values ****
#     Ncovar.epsilon - number of covariates for logit(epsilon)
#     covar.epsilon - matrix of covariates (Nsites x (Nyears-1) x N.covar.epsilon)

#     Ncovar.gamma - number of covariates for logit(gamma)
#     covar.gamma - matrix of covariates (Nsites x (Nyears-1) x N.covar.gamma)

#     Npred.psi  - number of predicted values needed for psi
#     pred.covar.psi - matrix of covariates where you want predictions to occur

#     Ncovar.p - number of covariates for logit(p)
#     covar.p - matrix of covariates (Nsites.visits x N.covar.p)

# 
cat(file="model.txt", "
    ############################################################
    
    model {
    # estimate the initial occupancy probabilities
    for(i in 1:Nsites){
       logit(psi[i,1]) = inprod(beta.psi[1:(Ncovar.psi-1)], covar.psi[i, 2:Ncovar.psi])
    }

    # estimate the extinction probabilities
    for(i in 1:(Nsites*(Nyears-1))){
      #epsilon[   site            , year  ]             <-               beta values                    covariates
       logit(epsilon[ covar.epsilon[i,1], covar.epsilon[i,2]]) <- inprod(beta.epsilon[1:(Ncovar.epsilon-2)], covar.epsilon[i,3:Ncovar.epsilon])
    }

    # estimate the gamma probabilities
    for(i in 1:(Nsites*(Nyears-1))){
      #gamma[      site            , year  ]     <-               beta values                    covariates
      logit(gamma[ covar.gamma[i,1], covar.gamma[i,2]]) <- inprod(beta.gamma[1:(Ncovar.gamma-2)], covar.gamma[i,3:Ncovar.gamma])
    }
    
    # project the psi forward for each site
    for(i in 1:Nsites){
       for(j in 2:Nyears){
         psi[i,j] <- psi[i,j-1]*(1-epsilon[i,j-1]) + (1-psi[i,j-1])*gamma[i,j-1]
       }
    }

    # set up the initial actual occupancy state model, i.e. is the site actually occupied or not
    for(i in 1:Nsites){
        z[i,1] ~  dbern(psi[i,1])
    }
    for(i in 1:Nsites){
       for(j in 2:Nyears){
         z[i,j] ~ dbern( z[i,j-1]*(1-epsilon[i,j-1]) + (1-z[i,j-1])*gamma[i,j-1])
       }
    }
    
    # the observation model.
    for(j in 1:Nsites.visits){
       logit(p.detect[j]) <- inprod(beta.p[1:(Ncovar.p-2)], covar.p[j, 3:Ncovar.p])
       p.z[j] <- z[Site[j],Year[j]]*p.detect[j]
       History[j] ~ dbern(p.z[j])
    }

    # priors
    #beta.psi[1] ~ dnorm(0, .01)  # this is the intercept
    for(i in 1:(Ncovar.psi-1)){
       beta.psi[i]   ~ dnorm(0, .0001)
    }

    #beta.epsilon[1] ~ dnorm(0, .01)  # this is the intercept
    for(i in 1:(Ncovar.epsilon-2)){
       beta.epsilon[i]   ~ dnorm(0, .0001)
    }

    #beta.gamma[1] ~ dnorm(0, .01) # this is the intercept
    for(i in 1:(Ncovar.gamma-2)){
       beta.gamma[i]   ~ dnorm(0, .0001)
    }

    #beta.p[1] ~ dnorm(0, .01) # this is the intercept on the logit scale
    for(i in 1:(Ncovar.p-2)){
       beta.p[i]   ~ dnorm(0, .0001)
    }
    
    # derived variables
    # number of occupied sites
    for(i in 1:Nyears){
       occ.sites[i] <- sum(z[1:Nsites,i])
    }

    # population growth rates
    for(i in 1:Nsites){
       lambda      [i,1:(Nyears-1)] <- psi[i,2:Nyears]/psi[i,1:(Nyears-1)]
       log.lambda.prime[i,1:(Nyears-1)] <- logit(psi[i,2:Nyears]) - logit(psi[i,1:(Nyears-1)]) 
       for(j in 1:(Nyears-1)){  # exp() cannot be vectorized
          lambda.prime[i,j] <- exp( log.lambda.prime[i,j] ) 
       }
    }

}
") # End of the model



# Next create the data.txt file.
# Initialize the data values using standard R code by either reading
# in from an external file, or plain assignment.

# get the data read in
# Data for detections should be a data frame with rows corresponding to sites
# and columns to visits.
# The usual 1=detected; 0=not detected; NA=not visited is used.
# Each year must have the same number of visit. Pad with missing if this is not true.

input.history <- read.csv(file.path("..","NSO.csv"), header=FALSE, skip=2, na.strings="-")
input.history$V1 <- NULL # drop the site number

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.history)
ncol(input.history)
range(input.history, na.rm=TRUE)
sum(is.na(input.history))

# Five years with 8 visits. input.history must have same number of visits/year
Nsites = nrow(input.history)
Nyears = 5
Nvisits= ncol(input.history)/Nyears


# Create the site covariates
site.covar <- data.frame(Site=1:Nsites, Browse=rep(c("Y","N"), length.out=Nsites))  # browse is a made up covariate
head(site.covar)

# Create the Yearly site (site x time) covariates. 
# yearlySiteCovs is a data frame with Nsites x Nyears rows which are in year-major, site-minor order.
# I.e. list covariates for year 1, sites 1... Nsites, then Site 2 for years .... etc
# The extinction and colonization rates will only depend on years 1... Nyears-1
yearsite.covar <- data.frame(Site   =rep(1:Nsites,  Nyears),
                             Year   =rep(1:Nyears,each=Nsites),
                             site.covar)
# We need to add the site covariates here as well
head(yearsite.covar)

# Create the Observation covarates (site x time x visit)
# Observation covaraites is a data frame with Nsite x Ntime x Nvisit rows
# which are ordered by columns in the design matrix 
# The order is IMPORTANT! Year, Visit, Site, i.e. stack the columns of matrix of visit covariates
obs.covar <- data.frame( Site   = rep(1:Nsites, Nyears * Nvisits),
                         Year   = rep(1:Nyears,  each=Nvisits*Nsites),
                         Visit  = rep(rep(1:Nvisits, each=Nsites), Nyears))
obs.covar <- merge(obs.covar, yearsite.covar, all.x=TRUE)  # add site-year covariates
obs.covar <- obs.covar[ order( obs.covar$Year, obs.covar$Visit, obs.covar$Site),]
obs.covar[1:30,]


History <- as.vector(unlist(input.history))  # stacks the columns
Site    <- obs.covar[,"Site"]
Year    <- obs.covar[,"Year"]
Visit   <- obs.covar[,"Visit"]


# Create design matrix for initial occupancy covariates
# ******* If you want year effects, use as.factor(Year) to get categorical variables rather than regression ******
covar.psi <- cbind(Site=site.covar[,"Site"], model.matrix( ~1, data=site.covar))
cbind(site.covar, covar.psi)[1:10,]

# Create design matrix for local extinction rate using yearsite covariates
covar.epsilon <- cbind(yearsite.covar[,c("Site","Year")], model.matrix( ~1, data=yearsite.covar))
cbind(yearsite.covar, covar.epsilon)[1:15,]

# Create design matrix for local colonization rate using yearsite covariates
covar.gamma <- cbind(yearsite.covar[,c("Site","Year")],model.matrix( ~1, data=yearsite.covar))
cbind(yearsite.covar, covar.gamma)[1:15,]

# Create design matrix for detection probabilities
covar.p <- cbind(Site, Year, model.matrix( ~as.factor(Year), data=obs.covar))
cbind(obs.covar, covar.p)[1:10,]

# Remove rows corresponding to missing values in all input data
no.visit <- is.na(History)

Site   <- Site   [ !no.visit]
Year   <- Year   [ !no.visit]
Visit  <- Visit  [ !no.visit]
History<- History[ !no.visit]
covar.p<- covar.p[ !no.visit,]
# we don't allow any missing values in the site, gamma, epsilon covariates.

Nsites.visits <- length(History) # after removal of missing values.

Ncovar.psi    <- ncol(covar.psi)
Ncovar.epsilon<- ncol(covar.epsilon)
Ncovar.gamma  <- ncol(covar.gamma)
Ncovar.p      <- ncol(covar.p)


# The datalist will be passed to JAGS with the names of the data
# values.
data.list <- list("Nsites","Nyears", "Nvisits", "Nsites.visits",
                  "History", "Site", "Year","Visit",
                  "Ncovar.psi",    "covar.psi",
                  "Ncovar.epsilon","covar.epsilon",
                  "Ncovar.gamma"  ,"covar.gamma",
                  "Ncovar.p"      ,"covar.p")

# check the list
data.list


# Next create the initial values.
# If you are using more than one chain, you need to create a function
# that returns initial values for each chain.

# We define the initial value of z as 1 if any visit resulted in a detection, other wise 0
init.z.long <- plyr::ddply(data.frame(History, Site, Year), c("Site","Year"), plyr::summarize, z=max(History))
init.z <- matrix(0, nrow=Nsites, ncol=Nyears)
init.z [ as.matrix(init.z.long[,c("Site","Year")])] <- init.z.long[,"z"]
init.z

# we will start at the same initial starting point for each chain even though this
# is not recommended. 
temp <- list(z=init.z, 
             beta.psi    =rep(0, Ncovar.psi-1),
             beta.epsilon=rep(0, Ncovar.epsilon-2),
             beta.gamma  =rep(0, Ncovar.gamma-2),
             beta.p      =rep(0, Ncovar.p-2)
)
init.list <- list(temp, temp, temp)



# Next create the list of parameters to monitor.
# The deviance is automatically monitored.
# 
monitor.list <- c("z",
                   "beta.psi", "beta.epsilon", "beta.gamma","beta.p",
                   "psi", "epsilon", "gamma", "p.detect", "occ.sites", 
                   "lambda", "lambda.prime") # parameters to monitor

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
str(results$BUGSoutput$summary)
results$BUGSoutput$summary[1:10,c("mean", "sd", "2.5%","97.5%","Rhat", "n.eff")]
results$BUGSoutput$summary[1:10,c("mean", "sd")]
results$BUGSoutput$summary[grepl("[1,",rownames(results$BUGSoutput$summary),fixed=TRUE),c("mean", "sd")]
results$BUGSoutput$summary[grepl("[2,",rownames(results$BUGSoutput$summary),fixed=TRUE),c("mean", "sd")]
results$BUGSoutput$summary[grepl("[1]",rownames(results$BUGSoutput$summary),fixed=TRUE),c("mean", "sd"),drop=FALSE]
results$BUGSoutput$summary[!grepl("[" ,rownames(results$BUGSoutput$summary),fixed=TRUE),c("mean", "sd"),drop=FALSE]
results$BUGSoutput$summary[grepl("lambda[1," ,rownames(results$BUGSoutput$summary),fixed=TRUE),c("mean", "sd")]
results$BUGSoutput$summary[grepl("lambda.prime[1," ,rownames(results$BUGSoutput$summary),fixed=TRUE),c("mean", "sd")]


# get just the means
results$BUGSoutput$mean
results$BUGSoutput$mean$psi

# the results$BUGSoutput$sims.array is a 3-d object [iterations, chains, variables]
dim(results$BUGSoutput$sims.array)
#results$BUGSoutput$sims.array[1:5,,]
results$BUGSoutput$sims.array[1:5,1,"psi[1,1]"]


# the results$BUGSoutput$sims.matrix is a 2-d object [iterations, variables] with chains stacked
# on top of each other
dim(results$BUGSoutput$sims.matrix)
results$BUGSoutput$sims.matrix[1:5,]
results$BUGSoutput$sims.matrix[1:5,"psi[1,1]"]


# make a posterior density plot
plotdata <- data.frame(parm=results$BUGSoutput$sims.matrix[,"psi[1,1]"], stringsAsFactors=FALSE)
head(plotdata)
postplot.parm <- ggplot2::ggplot( data=plotdata, aes(x=parm, y=..density..))+
  geom_histogram(alpha=0.3)+
  geom_density()+
  ggtitle("Posterior density plot for psi")
postplot.parm
ggsave(plot=postplot.parm, file='psi-posterior-p-t-psi-dot.png', h=4, w=6, units="in", dpi=300)



# make a trace plot (notice we use the sims.array here)
plotdata <- data.frame(psi=results$BUGSoutput$sims.array[,,"psi[1,1]"], stringsAsFactors=FALSE)
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
acf.parm <-acf( results$BUGSoutput$sims.matrix[,"psi[1,1]"], plot=FALSE)
acf.parm
acfplot.parm <- ggplot(data=with(acf.parm, data.frame(lag, acf)), aes(x = lag, y = acf)) +
  ggtitle("Autocorrelation plot for psi")+
  geom_hline(aes(yintercept = 0)) +
  geom_segment(aes(xend = lag, yend = 0))
acfplot.parm
ggsave(plot=acfplot.parm, file="psi-acf-p-t-psi-dot.png",h=4, w=6, units="in", dpi=300)


# Model assessment (Bayesian Posterior Predictive Checks) not yet implemented here.
