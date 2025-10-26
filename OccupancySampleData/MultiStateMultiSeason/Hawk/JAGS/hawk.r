# Multi-state occupancy model using the psi0, r0, psi, r, p and deltat parameterization


# phi0  = initial occupancy probabilities in each state parameterized as  phi0[1]*(1-phi0[1]), phi0[1](phi0[2]])

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



# Ensure that the R2jags package is available.

library("R2jags")  # used for call to JAGS
library(coda)
library(ggplot2)
library(reshape2)

# The model file.
# The cat() command is used to save the model to the working directory.
# Notice that you CANNOT have any " (double quotes) in the bugs code
# between the start and end of the cat("...",) command.

cat(file="model.txt", "
############################################################
var  psi0.dm [nsites, psi0.nbeta], psi0.beta [psi0.nbeta],
     r0.dm   [nsites, r0.nbeta],   r0.beta   [r0.nbeta],
     Cpsi.dm [nsites*(nyears-1)*3, Cpsi.nbeta], Cpsi.beta[Cpsi.nbeta],
     Cr.dm   [nsites*(nyears-1)*3, Cr.nbeta],   Cr.beta[Cr.nbeta],
     site.year.state[nsites*(nyears-1)*3, 3], 
     TR      [nsites, nyears-1, 3, 3],        # conditional transition matrices
     p1.dm   [nobsocc, p1.nbeta],  p1.beta   [p1.nbeta],
     p2.dm   [nobsocc, p2.nbeta],  p2.beta   [p2.nbeta],
     delta.dm[nobsocc,delta.nbeta],delta.beta[delta.nbeta],
     P      [nobsocc, 3, 3],
     Occ    [nyears, nsites];

model {
   # Input data consists of

   #   nsites    - number of sites in the occupancy model
   #   nyears    - number of years in the study

   #   ostate    - observed occupancy state (incremented by 1 because JAGS indexes starting at 1 rather than 0
   #                ostate=1 = no observed occupancy; =2 is observed in state 1; =3 is observed in state 3)
   #   ostate.site - site number of observed occupancy state
   #   ostate.year - year of observed occupancy state
   #   ostate.survey - survey (within a year) of observed occupancy state
   #   nobsocc       - number of observed occupancy states

   #   psi0.dm   - design matrix for logit(psi(0)) (initial occupancy probability (state 1 or state 2))
   #   psi0.nbeta- number of columns in psi0 dm (i.e # of parameters)

   #   r0.dm     - design matrix for logit(r(0))   (initial probability of being in state 2 | in state 1)
   #   r0.nbeta  - number of columns in r0 dm   (i.e. # of parameters)

   #   Cpsi.dm    - design matrix for conditional probability of occupancy in year t+1 | current state in year t
   #   Cpsi.nbeta - number of columns in Cpsi.dm

   #   Cr.dm      - design matrix for conditional probability of breeding in year t+1 | current state in year t
   #   Cr.nbeta   - number of columns in Cr.dm

   #   p1.dm     - design matrix for logit(p1(obsocc)) (probability of detection | actual state is 1 )
   #   p1.nbeta  - number of columns in p1 dm
 
   #   p2.dm     - design matrix for logit(p2(obsocc)) (probability of detection | actual state is 2 )
   #   p2.nbeta  - number of columns in p2 dm
 
   #   delta.dm  - design matrix for logit(delta(obsocc)) (probability of detect in state 2 | actual state is 2)
   #   delta.nbeta- number of columns in delta.dm

   
   #  compute the initial occupancy and r for each site
   for(i in 1:nsites){
      logit(psi0[i]) <- inprod(psi0.dm[i,], psi0.beta[])
      logit(r0  [i]) <- inprod(r0.dm[i,]  , r0.beta[])
   }

   #  compute the p1, p2, and delta
   for(i in 1:nobsocc){
      logit(p1   [i])    <-  inprod(p1.dm[i,],     p1.beta[])
      logit(p2   [i])    <-  inprod(p2.dm[i,],     p2.beta[])
      logit(delta[i])    <-  inprod(delta.dm[i,],  delta.beta[])

      # now for P probability of detection in state i if true state is j
      # index for P[ observed occupancy, true state, obs state]
      #                                                    Observed State                                        
      #                                    0                    1                           2 
      #------------------------------------------------------------------------------------------------
      #True State           0|             1                    0                           0         
      #                     1|             1 - p[1]            p[1]                         0         
      #                     2|             1 - p[2]            p[2]*(1 - Delta[2,2])    p[2]*Delta[2,2]
      P[i, 0+1, 0+1] <- 1  # don't forget that JAGS index from 1 rather than 9
      P[i, 0+1, 1+1] <- 0
      P[i, 0+1, 2+1] <- 0
      P[i, 1+1, 0+1] <- 1-p1[i]
      P[i, 1+1, 1+1] <- p1[i]
      P[i, 1+1, 2+1] <- 0
      P[i, 2+1, 0+1] <- 1-p2[i]
      P[i, 2+1, 1+1] <- p2[i]*(1-delta[i])
      P[i, 2+1, 2+1] <- p2[i]*delta[i]
   }

   #  compute phi0 (for each site) which is p(occupancy in state 0, 1, or 2)
   phi0[0+1, 1:nsites] = 1 - psi0[1:nsites]                # state 0. Don't forget that JAGS starts at index 1
   phi0[1+1, 1:nsites] = psi0[1:nsites] * (1-r0[1:nsites]) # state 1
   phi0[2+1, 1:nsites] = psi0[1:nsites] * r0[1:nsites]     # state 2

   #  compute the transition matrix (site, year,)
   # The parameters Psi[0], Psi[1], Psi[2], R[0,2], R[1,2], and R[2,2] are used to construct 
   # the Phi_t matrix for transitions each interval (midway down second column page 826 of MacKenzie et al. 2009):
   #  
   #  1 - Psi[0](t)          Psi[0](t)*{1 - R[0,2](t)}          Psi[0](t)*R[0,2](t)
   #  1 - Psi[1](t)          Psi[1](t)*{1 - R[1.2](t)}          Psi[1](t)*R[1,2](t)
   #  1 - Psi[2](t)          Psi[2](t)*{1 - R[2,2](t)}          Psi[2](t)*R[2,2](t)
   for(i in 1:(nsites*(nyears-1)*3)){
      logit(Cpsi[i]) = inprod(Cpsi.dm[i,],     Cpsi.beta[]  )
      logit(Cr[i]  ) = inprod(Cr.dm  [i,],       Cr.beta[]  )

      TR[ site.year.state[i,1], site.year.state[i,2], site.year.state[i,3], 0+1] <- 1- Cpsi[i]
      TR[ site.year.state[i,1], site.year.state[i,2], site.year.state[i,3], 1+1] <- Cpsi[i]*(1-Cr[i])
      TR[ site.year.state[i,1], site.year.state[i,2], site.year.state[i,3], 2+1] <- Cpsi[i]*Cr[i]
   }


   # model initial occupancy latent state  Occ[year, site]
   for (i in 1:nsites) {
	  	Occ[1,i] ~ dcat(phi0[(0:2)+1, i]) 
   }

   # model subsequent occupancy latent states. 
   for(i in 1:nsites){
       for(j in 2:nyears){
          Occ[j,i] ~dcat( TR[i, j-1, Occ[j-1,i], (0:2)+1])
          #Occ[j,i]  ~dcat( c(.25, .50, .25))

       }
   }

   # now for the observed data. 
   for(i in 1:nobsocc){
       ostate[i] ~ dcat( P[i, Occ[ostate.year[i], ostate.site[i]], 1:3])
   }



   # priors for the betas
   for(i in 1:psi0.nbeta){
      psi0.beta[i] ~ dnorm(0, .0001)  # mostly flat prior
   }
   for(i in 1:r0.nbeta){
      r0.beta[i]   ~ dnorm(0, .0001)  # mostly flat prior
   }
   for(i in 1:Cpsi.nbeta){
      Cpsi.beta[i]   ~ dnorm(0, .0001)  # mostly flat prior
   }
   for(i in 1:Cr.nbeta){
      Cr.beta[i]   ~ dnorm(0, .0001)  # mostly flat prior
  }
   for(i in 1:p1.nbeta){
      p1.beta[i]   ~ dnorm(0, .0001)  # mostly flat prior
   }
   for(i in 1:p2.nbeta){
      p2.beta[i]   ~ dnorm(0, .0001)  # mostly flat prior
   }
   for(i in 1:delta.nbeta){
      delta.beta[i]   ~ dnorm(0, .0001)  # mostly flat prior
   }
}
") # End of the model

#



# Next create the data.txt file.
# Initialize the data values using standard R code by either reading
# in from an external file, or plain assignment.

input.data <- read.csv(file.path("..", "hawk.csv"), header=TRUE, as.is=TRUE, strip.white=TRUE)

detect.history <- input.data[, paste("V",1:54,sep="")]
unit.cov       <- input.data[, "tsh", drop=FALSE]
# no survey level covaraites

unit.cov$site <- 1:nrow(unit.cov)

nsites  <- nrow(detect.history)
nyears  <- 6
nsurveys<- 9

# create the detection history matrix - need to include the site, year, survey, and state value
# Note that we increase the state by 1 because JAGS indexes from 1 rather than from 0
# so ostate=1 implies observed stated = 0; ostate=2 implies observed state = 1, etc
obsocc     <- data.frame(ostate  = as.vector(as.matrix(detect.history+1)),
                         site    = rep(1:nsites,   length.out=nsites*nyears*nsurveys),
                         year    = rep(1:nyears,   each=nsites*nsurveys),
                         survey   =rep(rep(1:nsurveys, each=nsites), length.out=nsites*nyears*nsurveys),
                         stringsAsFactors=FALSE)
head(obsocc)
# remove all missing values
select <- !is.na(obsocc$ostate)
obsocc <- obsocc [select,]

nobsocc <- nrow(obsocc)

# create the design matrix for the initial psi0 and r0 vectors.
psi0.dm <- model.matrix( ~tsh, data=unit.cov )
psi0.nbeta <- ncol(psi0.dm)

r0.dm <- model.matrix( ~1, data=unit.cov )
r0.nbeta <- ncol(r0.dm)

# Create the design matrix for the Conditional psi and r .
# We need to construct a pseudo-observed occupancy with 3 rows for every actual row corresponding to
# state 0, 1, and 2 to allow the psi and r terms to vary with state.
site.year.state <- data.frame( site   =rep(1:nsites, each=3, length.out=nsites*(nyears-1)*3),
                               year   =rep(1:(nyears-1), each=nsites*3),
                               state  =rep((0:2)+1,  length.out=nsites*(nyears-1)*3))
Cpsi.dm    <- model.matrix(~as.factor(state), data=site.year.state)
Cpsi.nbeta <- ncol(Cpsi.dm)

Cr.dm      <- model.matrix(~as.factor(state), data=site.year.state)
Cr.nbeta   <- ncol(Cr.dm)

# create the design matrix for p1, p2, and delta
p1.dm <- model.matrix(~1, data=obsocc)
p1.nbeta <- ncol(p1.dm)

p2.dm <- model.matrix(~1, data=obsocc)
p2.nbeta <- ncol(p2.dm)

delta.dm <- model.matrix(~1, data=obsocc)
delta.nbeta <- ncol(delta.dm)


# The datalist will be passed to JAGS with the names of the data values


data.list <- list(nsites=nsites, nyears=nyears,
                  ostate=obsocc$ostate, ostate.site=obsocc$site, 
                                 ostate.year=obsocc$year, 
                                 ostate.survey=obsocc$survey, nobsocc=nobsocc,  # observed occupancy values
                  psi0.dm=psi0.dm, psi0.nbeta=psi0.nbeta,
                  r0.dm=r0.dm,     r0.nbeta  =r0.nbeta,
                  Cpsi.dm=Cpsi.dm, Cpsi.nbeta=Cpsi.nbeta,
                  Cr.dm=Cr.dm,     Cr.nbeta=Cr.nbeta, site.year.state=as.matrix(site.year.state),
                  p1.dm  =p1.dm,   p1.nbeta=p1.nbeta,
                  p2.dm  =p2.dm,   p2.nbeta=p2.nbeta,
                  delta.dm=delta.dm,delta.nbeta=delta.nbeta)

# check the list
#data.list


# Next create the initial values.

# Initialize the actual state as the observed occupancy

init.occ <- plyr::ddply(obsocc, c("site","year"), plyr::summarize,
                   iocc=max(ostate))
head(init.occ)
Occ <- matrix(3, nrow=nyears, ncol=nsites)
#Occ[ cbind(init.occ$year, init.occ$site)] <- init.occ$iocc
max(Occ, na.rm=TRUE)

init.list <- list(
      list(Occ=Occ),
      list(Occ=Occ),
      list(Occ=Occ)
      )  # end of list of lists of initial value



# Next create the list of parameters to monitor.
# The deviance is automatically monitored.
# 
monitor.list <- c("psi0.beta","psi0",
                  "r0.beta","r0", 
                  "Cpsi.beta", "Cpsi",
                  "Cr.beta","Cr",
                  "p1.beta", "p1",
                  "p2.beta", "p2",
                  "delta.beta","delta",
                  "phi0","P","TR",
                  "Occ")        # latent true state (incremented by 1)
 

   


# Finally, the actual call to JAGS
set.seed(234234)  # intitalize seed for MCMC 

results <- jags( 
      data      =data.list,   # list of data variables
      inits     =init.list,   # list/function for initial values
      parameters=monitor.list,# list of parameters to monitor
      model.file="model.txt",  # file with bugs model
      n.chains=3,
      n.iter  =8000,          # total iterations INCLUDING burn in
      n.burnin=2000,          # number of burning iterations
      n.thin=2,               # how much to thin
      DIC=TRUE,               # is DIC to be computed?
      working.dir=getwd()    # store results in current working directory
      )


save(list=c("results"), file="results.Rdata")
#load(file="results.Rdata")

extract.index<- function(names){
    # extract the indices from the names, i.e. if names is psi0[10] return 10, if T[1,3,2] return 1, 2, 3
    # get rid of everything up to first [
    names <- substring(names, regexpr("[", names, fixed=TRUE)+1)
    # get rid of last ]
    names <- gsub("]","", names, fixed=TRUE)
    res <- strsplit(names, ",", fixed=TRUE)
    res <- plyr::laply(res, function (x){
        as.numeric(x)
    })
    res
}



# get the summary table
results$BUGSoutput$summary
results$BUGSoutput$summary[,c("mean", "sd", "2.5%","97.5%","Rhat", "n.eff")]

##############################################################
# Get the initial probability of occupancy.
select <- grepl("psi0[",row.names(results$BUGSoutput$summary), fixed=TRUE)
sum(select)
row.names(results$BUGSoutput$summary)[select]
psi0.res <- as.data.frame(results$BUGSoutput$summary[select,], stringsAsFactors=FALSE)
psi0.res$lcl <- psi0.res$"2.5%"
psi0.res$ucl <- psi0.res$"97.5%"
psi0.res$site <- extract.index(row.names(results$BUGSoutput$summary)[select])
psi0.res <- merge(psi0.res, unit.cov, by="site")
psi0.res[1:10,] 
  
ggplot(data=psi0.res, aes(x=tsh, y=mean))+
   ggtitle("Initial occupancy probability")+
   geom_point()+
   geom_ribbon(aes(ymin=lcl, ymax=ucl), alpha=0.2)+
   ylim(0,1)+ylab("Initial occupancy probability")+xlab("Area of suitable habitat ")

select <- grepl("r0[",row.names(results$BUGSoutput$summary), fixed=TRUE)
sum(select)
row.names(results$BUGSoutput$summary)[select]
results$BUGSoutput$summary[select,][1:10,]


##############################################################
# Get the initial three initial states.
select <- grepl("phi0[",row.names(results$BUGSoutput$summary), fixed=TRUE)
sum(select)
row.names(results$BUGSoutput$summary)[select]
phi0.res <- as.data.frame(results$BUGSoutput$summary[select,], stringsAsFactors=FALSE)
phi0.res$lcl <- phi0.res$"2.5%"
phi0.res$ucl <- phi0.res$"97.5%"
phi0.res$state <-  extract.index(row.names(results$BUGSoutput$summary)[select])[,1]-1  # don't forget that JAGS had to add 1 to state variable
phi0.res$site  <-  extract.index(row.names(results$BUGSoutput$summary)[select])[,2]  
phi0.res <- merge(phi0.res, unit.cov, by="site")
phi0.res[1:10,] 

ggplot(data=phi0.res, aes(x=tsh, y=mean))+
  ggtitle("Initial probability of being in states")+
  geom_point()+
  geom_ribbon(aes(ymin=lcl, ymax=ucl), alpha=0.2)+
  ylim(0,1)+ylab("Probability")+xlab("Area of suitable habitat at breeding scale")+
  facet_wrap(~state, ncol=2)

##############################################################
# Get the estimated transitions matrix
select <- grepl("TR",row.names(results$BUGSoutput$summary), fixed=TRUE)
sum(select)
row.names(results$BUGSoutput$summary)[select]

select <- grepl("TR[2,1",row.names(results$BUGSoutput$summary), fixed=TRUE)
sum(select)
extract.index(row.names(results$BUGSoutput$summary)[select])

round(matrix(results$BUGSoutput$summary[select,"mean"], nrow=3, ncol=3),3)
round(matrix(results$BUGSoutput$summary[select,"sd"  ], nrow=3, ncol=3),3)

##############################################################
# Estimates of catchability
select <- grepl("p1[",row.names(results$BUGSoutput$summary), fixed=TRUE)
sum(select)
cbind(obsocc,results$BUGSoutput$summary[select,])[1:10,]

select <- grepl("p2[",row.names(results$BUGSoutput$summary), fixed=TRUE)
sum(select)
cbind(obsocc,results$BUGSoutput$summary[select,])[1:10,]

select <- grepl("delta[",row.names(results$BUGSoutput$summary), fixed=TRUE)
sum(select)
cbind(obsocc,results$BUGSoutput$summary[select,])[1:10,]


##############################################################
# Estimates of posterior probability of occupancy in each state
# We need to use the Occ[] array of states and them compute various averages
dim(results$BUGSoutput$sims.matrix)
select <- grepl("Occ[",colnames(results$BUGSoutput$sims.matrix), fixed=TRUE)
sum(select)

occ <- reshape2::melt(as.data.frame(results$BUGSoutput$sims.matrix[,select]-1),
                    variable.name="year.site",
                    value.name   ="Occupancy")
temp <- extract.index(occ$year.site)
occ$site <- temp[,2]
occ$year <- temp[,1]
occ$Occ1 <- occ$Occupancy %in% c(1,2)  # level 1 occupancy
occ$Occ2 <- occ$Occupancy==2  # level 2 occupancy
rm(temp)
head(occ)

mean.occ.site <- plyr::ddply(occ, c("year","site"), plyr::summarize, 
                        Occ1 = mean(Occ1),
                        Occ2 = mean(Occ2))
mean.occ.site <- reshape2::melt(mean.occ.site,
                                id.var=c('year',"site"),
                                variable.name='OccupancyLevel',
                                value.name="prob")
mean.occ.site


head(mean.occ.site)
mean.occ <- plyr::ddply(mean.occ.site, c("year","OccupancyLevel"), plyr::summarize, 
                             prob.mean=mean(prob),
                             prob.sd  =sd(prob),
                             prob.n   =length(prob))

mean.occ$lcl <- mean.occ$prob.mean - 2*mean.occ$prob.sd/sqrt(mean.occ$prob.n)
mean.occ$ucl <- mean.occ$prob.mean + 2*mean.occ$prob.sd/sqrt(mean.occ$prob.n)

ggplot(data=mean.occ, aes(x=year, y=prob.mean, color=OccupancyLevel))+
  ggtitle("Estimated occupancy levels over time")+
  geom_line( position=position_dodge(width=.1))+
  geom_errorbar( aes(ymin=lcl, ymax=ucl), width=.2,  position=position_dodge(width=.1))+
  ylim(0,1)+ylab("Estimated probability")+
  xlab("Year")
