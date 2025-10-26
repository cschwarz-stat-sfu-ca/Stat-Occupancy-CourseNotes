# Single Species Multi Method Single Season Occupancy
# Coded 2024-03-20 Carl James Schwarz


# Other notes
#  (1) in the AIC table, use "use.aicc=TRUE' to get AICc rather than the simpler AIC.
#  (2) the occMod_mm() function uses three covariates data frames
#       psi.cov  for covariates on sites
#       theta.cov for covariates on the theta parameters (site x visit).
#          Unfortunately, you need to define a "VISIT" covariate here because RPresence does not
#       p.cov     for covariates on the p parameter (site x visit x method)
#          RPresence defined SURVEY and DEVICE for these covariates as needed
#  (3) modAvg doesn't work for psi for some reason. use my.modAvg() function until JH fixes

my.modAvg <- function (aic.tab, param = "psi", parmgrp = "real", index = 1:nrow(aic.tab$table), 
          conf = 0.95, predict = FALSE, newdata, trace=FALSE,replaceNaN.with.zero=TRUE) 
{
  linkfn = "logit"
  if (length(grep("phi|nu", param)) > 0) 
    linkfn = "exp"
  aic.table <- list(table = aic.tab$table[index, ], models = aic.tab$models[index])
  class(aic.table) <- "aic.table"
  jaic = grep("AIC", colnames(aic.table$table))[1]
  jdaic = grep("AIC", colnames(aic.table$table))[2]
  aic.table$table$DAIC <- aic.table$table[, jaic] - min(aic.table$table[, 
                                                                        jaic])
  aic.table$table$modlike <- exp(-aic.table$table[, jdaic]/2)
  aic.table$table$wgt <- aic.table$table$modlike/sum(aic.table$table$modlike)
  if (!predict) {
    i = which(names(aic.table$models[[1]]) == parmgrp)
    nr = nrow(get(param, aic.table$models[[1]][[i]]))
    est <- array(unlist(lapply(aic.table$models, function(xx) get(param, 
                                                                  xx[[i]])[, 1])), dim = c(nr, length(aic.table$models)))
    se <- array(unlist(lapply(aic.table$models, function(xx) get(param, 
                                                                 xx[[i]])[, 2])), dim = c(nr, length(aic.table$models)))
  
  }
  else {
    pred <- lapply(aic.table$models, function(xx) predict(xx, 
                                                          newdata, param = param))
    est <- array(unlist(lapply(pred, function(xx) xx[, 1])), 
                 dim = c(nrow(newdata), length(aic.table$models)))
    se <- array(unlist(lapply(pred, function(xx) xx[, 2])), 
                dim = c(nrow(newdata), length(aic.table$models)))
  }
  if(replaceNaN.with.zero){ # CJS replace NaN in se with 0. Typically occur if estimate is 0 or 1
     se[is.na(se)] <- 0
  }
  ma <- est %*% aic.table$table$wgt
  ma.se <- sqrt((se^2 + (est - as.vector(ma))^2) %*% aic.table$table$wgt)
  if (linkfn == "logit") {
    logit.est <- log(ma/(1 - ma))
    logit.se <- ma.se
    on.bounds <- ma %in% 1 | ma %in% 0
    logit.se[on.bounds] <- 0
    logit.se[!on.bounds] <- ma.se[!on.bounds]/(ma[!on.bounds] * 
                                                 (1 - ma[!on.bounds]))
  }
  else {
    logit.est = log(ma)
    logit.se = ma.se
    logit.se[ma == 0] = 0
  }
  alpha <- (1 - conf)/2
  z <- -qnorm(alpha)
  lower <- sapply(z, function(zz) logit.est - zz * logit.se)
  #browser()
  lower <- as.matrix(lower)  # CJS avoid problems with a single parameter value
  colnames(lower) <- paste("lower", conf, sep = "_")
  upper <- sapply(z, function(zz) logit.est + zz * logit.se)
  upper <- as.matrix(upper)  # CJS ditto avoid problems with single parameter value
  colnames(upper) <- paste("upper", conf, sep = "_")
  if (linkfn == "logit") {
    lower <- plogis(lower)
    upper <- plogis(upper)
  }
  else {
    lower = exp(lower)
    upper = exp(upper)
  }
  result <- data.frame(est = ma, se = ma.se, lower, upper)
  if(trace)browser()
  if (!predict) {
    rownames(result) <- rownames(get(param, aic.table$models[[1]][[i]]))
  }
  else {
    rownames(result) <- rownames(newdata)
  }
  return(result)
}



#-------------------------------Setup pckgs and GE--------------------------------------
library(plyr)
library(readxl)
library(RPresence)
library(ggplot2)

# clear global environment
rm(list = setdiff(ls(), "my.modAvg"))

#-------------------------------Data Import & Checks------------------------------------
# get the data read in
# MultiMethod data structure for the capture history is
#  v1m1 v1m2 v2m1 v2m2 v3m1 v3m2
# where v1=visit 1, m1=method 1 etc
# So, for each visit, you have set of 0/1/. indicating if the method did not detect, did detect, was not used.
# Normally, each method is used on each visit. This does NOT have to happen every time, but there must be some visits where 
# multiple methods are used.


files <- dir(file.path(".."))
load(file.path("..", files[grep("rdata",files)]))
names(mmRMark)[[2]]
input.history <- mmRMark[[2]]
input.history <- strsplit(input.history$ch, split="")
input.history <- matrix(as.numeric(unlist(input.history)), ncol = 9, byrow = TRUE)

n.methods=3  # number of methods to be used
n.visits =ncol(input.history)/n.methods
n.sites  =nrow(input.history)

cat("Methods :", n.methods,";   visits =", n.visits,  "\n")
# check number of sites; number of visits)
nrow(input.history)
ncol(input.history) # this will be number of visits x number of methods 
range(input.history, na.rm=TRUE) # check that all values are either 0 or 1
sum(is.na(input.history)) # are there any missing values? This would indicate a method was not used on a visit

head(input.history)

site_covar <- data.frame(site.cov=rnorm(nrow(input.history)))  # no site level covariates

# summarize the number of records in each category from the newly created "SiteCat" field
#xtabs(~SiteCat, data=site_covar,exclude=NULL, na.action=na.pass)

head(site_covar)
str(site_covar)


#------------------------------Data Manipulation, Formatting-----------------------------
# p covariate. Here you stack the methods within visits across sites

# In this case the p.cov data.frame is constructed by stacking the covariate values

p.cov <- data.frame(site=as.factor(rep(1:n.sites , n.visits*n.methods)),
                         # this first bit above stacks the site field based on the input.history dataframe
                    VISIT =as.factor(rep(1:n.visits, each=n.sites*n.methods)),
                    device=as.factor(rep(rep(1:n.methods, each=n.sites), n.visits))
                 )
head(p.cov, n=25)

# Theta - probability of species being available if present in the site, varies by site and VISIT
# RPresence does not define the visit and so you must do it here.
theta.cov <- data.frame(site=as.factor(rep(1:n.sites, n.visits)),
                        VISIT=as.factor(rep(1:n.visits, each=n.sites)))
head(theta.cov, n=20)





#----------------------------Create the *.pao file for PRESENCE--------------------------
birds.pao <- RPresence::createPao(input.history,  # easier to define the covariate values in the call to occMod
                                 methods=n.methods,  
                                 title='birds')
birds.pao

# check the structure of the pao data - previously the RainIndex (a continuous covariate) had been converted to a factor
str(birds.pao)




#------------------------------------ Model Fit Approach 1 ---------------------------------------------------
# Parameters for Multi Method Single Season occupancy are
#.  psi - occupancy probability
#.  p  - probability of detection on visit (SURVEY) and method (DEVICE)
#.  theta - probability that individuals are available for detection using method x, given presence
#           theta is really used for multi-scale studies where local occupancy could be different.
#.          for multi-method, set theta to 1

mod.pdot_base <- RPresence::occMod(
  model=list(psi~1, p~1, theta~1),
  cov.list=list(psi.cov=site_covar, theta.cov=theta.cov, p.cov=p.cov),
  type="so.mm", 
  data=birds.pao)

summary(mod.pdot_base)

# get estimates on the real scale
mod.pdot_base$real$psi
mod.pdot_base$real$theta
mod.pdot_base$real$p


# returning the results on the standard/regular scale (converted back from logit)
fitted(mod.pdot_base, param="psi")
fitted(mod.pdot_base, param="p")
fitted(mod.pdot_base, param="theta")

names(mod.pdot_base)
mod.pdot_base$data
mod.pdot_base$dmat




# ------------------------------- Model Fit Approach 2 ---------------------------------------------------
# Newer code to define a fill list of models to fit
# Notice the commas between the column and the placement of the quotes
model.list.csv <- textConnection("
p,               psi,       theta
~1,              ~1,        ~1
~SURVEY,         ~1,        ~1
~1,              ~1,        ~VISIT
~DEVICE,         ~1,        ~1
~DEVICE,         ~1,        ~VISIT
~SURVEY,         ~1,        ~VISIT")
 


#~RainIndex, ~1
#~RainIndex, ~SiteCat
#")


model.list <- read.csv(model.list.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)
model.list


# fit the models
model.fits <- plyr::alply(model.list, 1, function(x,detect.pao){
  cat("\n\n***** Starting ", unlist(x), "\n")
  fit <- RPresence::occMod(model=list(as.formula(paste("psi",x$psi)),
                                      as.formula(paste("p"  ,x$p  )),
                                      as.formula(paste("theta", x$theta))),
                           cov.list=list(psi.cov  =site_covar,
                                         theta.cov=theta.cov,
                                         p.cov    =p.cov),  # notice we add covariate list for theta
                           data=detect.pao,type="so.mm")
  fit
},detect.pao=birds.pao)



# Look the output from a specific model
check.model <- 2

names(model.fits[[check.model]])
model.fits[[check.model]]$beta

names(model.fits[[check.model]]$real)
model.fits[[check.model]]$real$psi
model.fits[[check.model]]$real$p[1:9,]
model.fits[[check.model]]$real$theta

names(model.fits[[check.model]]$derived)
model.fits[[check.model]]$derived$psi_c[1:10,]
tail(model.fits[[check.model]]$derived$psi_c)
model.fits[[check.model]]$derived$p_c[1:10,]

#-------------------------------------Model averaging----------------------------------
AICc.table <- RPresence::createAicTable(model.fits, use.aicc=TRUE) # note to use AICc
AICc.table$table

names(AICc.table)

ma.p <- RPresence::modAvg(AICc.table, param="p")
ma.p[grepl("unit1$", rownames(ma.p)),]
my.modAvg(AICc.table, param="p", replaceNaN.with.zero=TRUE)

# ma.psi <- RPresence::modAvg(AICc.table, param="psi") # does not work
ma.psi <- my.modAvg(AICc.table, param="psi") # use instead
ma.psi

ma.theta <- RPresence::modAvg(AICc.table, param="theta")
ma.theta
ma.theta <- my.modAvg(AICc.table, param="theta", replaceNaN.with.zero=TRUE)
ma.theta

