# RMark additional functions

# 2018-02-15 Software contributed by Carl Schwarz (cschwarz.stat.sfu.ca@gmail.com)

RMark.random.occup.stack <- function(input.history, surveys.per.season){
   # Stack the input.history to get it ready for a random occupancy model
   # This takes a multi-season file and converts it to a series of individual years
   # For example, the history
   #      010 100
   # is converted to
   #      010 ...
   #      ... 100
   # so that separate psi can be fit for each season, but common p' across seasons can be fit.
   #
   require(plyr)
  
   n.season <- length(surveys.per.season)
   
   input.history <- input.history[rep(1:nrow(input.history), each=n.season),]
   input.history$Season <- 1:n.season  
    
   # now to "blank out" (set to missing) parts of the history outside of each season
   input.history <- plyr::ddply(input.history, "Season", function(x, surveys.per.season) {
     start <- cumsum(c(1,surveys.per.season))[x$Season[1]]
     end   <- cumsum(c(surveys.per.season,1))[x$Season[1]]
     extract <- substr(x$ch, start, end)
      x$ch    <- paste(rep(".",sum(surveys.per.season)), sep="", collapse="")
     substr(x$ch, start, end) <- extract
     x
   }, surveys.per.season=surveys.per.season)
}


##################################################################################
# 2019-12-23 CJS Need to compute the se for the derived psi as this is no longer computed

RMark.add.derived <- function(model){
  # Add additional derived variables to RMark models
  # This is particularly true for occupancy models where the time specific psi's can be in different places
  
  
  # RDOccupEG. RMark stores all of the psi's in the derived parameter 
  if("RDOccupEG" %in% class(model)){ # single species multi season psi, p, gamma, epsilon parameterization
    model$results$derived$all_psi <- model$results$derived$"psi Probability Occupied"
    model$results$derived$all_psi$se <- sqrt(diag(model$results$derived.vcv$"psi Probability Occupied"))
  }
  if("RDOccupPG" %in% class(model)){ # single species multi season psi, p, gamma,  parameterization
    model$results$derived$all_psi <- get.real(model, "Psi", se=TRUE)
  }
  if("RDOccupPE" %in% class(model)){ # single species multi season psi, p, epsilon,  parameterization
    model$results$derived$all_psi <- get.real(model, "Psi", se=TRUE)
  }
  
  
  model
}
  


##################################################################################
# Model averaging derived parameters
# http://www.phidot.org/forum/viewtopic.php?f=21&t=2865&p=9588&hilit=model+averaging+derived+estimates#p9588
# 2019-12-23 CJS add .drop=FALSE if model average only a single models

RMark.model.average.derived <- function(model.set, parameter, estimate="estimate", se="se", alpha=0.25){
  # model average derived parameters
  # like RMark, alpha=.025 gives a 95% confidence interval
  # The aic table has all of the model information, but also the aic table at the end
  
  require(plyr)
  require(boot)
  #browser()
  # Extract the estimates and standard errors
  est    <- plyr::laply(model.set[-length(model.set)], function(x){x$results$derived[[parameter]][,estimate]},.drop=FALSE)
  est.se <- plyr::laply(model.set[-length(model.set)], function(x){x$results$derived[[parameter]][,se      ]},.drop=FALSE)
  AICc   <- plyr::laply(model.set[-length(model.set)], function(x){x$results$AICc}) 
 
  # Get the model averaged values
  ma.est <- data.frame(model.average(list(estimate=est, se=est.se, AICc=AICc)))
  wald.params <- c("lambda Rate of Change","log odds lambda")
  if (!(parameter %in% wald.params)){
    ma.est$est.logit <- logit(ma.est$estimate)
    ma.est$est.logit.se <- ma.est$se/(ma.est$estimate*(1-ma.est$estimate))
    ma.est$ci.logit.lower <- ma.est$est.logit - qnorm(.975)*ma.est$est.logit.se
    ma.est$ci.logit.upper <- ma.est$est.logit + qnorm(.975)*ma.est$est.logit.se
    ma.est$lcl <- inv.logit(ma.est$ci.logit.lower)
    ma.est$ucl <- inv.logit(ma.est$ci.logit.upper)
    keep.params <- c("estimate","se","lcl","ucl")
    ma.est <- ma.est[,which(names(ma.est) %in% keep.params)]
  } 
  else {
    ma.est$lcl <- ma.est$estimate - qnorm(.975)*ma.est$se
    ma.est$ucl <- ma.est$estimate + qnorm(.975)*ma.est$se 
  }
  ma.est
  
}
