# Additional Functions for RPresence

# 2018-08-15 Code contributed by Carl James Schwarz (cschwarz.stat.sfu.cs@gmail.com)

RPresence.add.derived <- function(model){
  # Add additional derived variables to RPresence models
   
  # DO.1. RPresence stores the initial psi in the real parameters and then
  #       puts the psi in subsequent years in the derived parameters.
  #       The "all_psi" puts both together in the derived parameters to make it easier to work with the full set of psi's
  
  require(plyr)
  
  if("do1" %in% class(model)){ # single species multt season psi, p, gamma, epsilon parameterization
      # Combine the initial psi and subsequent psi into one derived parameters
      model$derived$all_psi <- rbind( model$real$psi, model$derived$psi)
      
      # Estimate the occupancy growth rate (lambda) for each site and year
      all_psi <- model$derived$all_psi
      all_psi$unit <- substr(row.names(all_psi),1,-1+regexpr("_",row.names(all_psi), fixed=TRUE))
      all_psi$time <- as.numeric(substring(row.names(all_psi),1+regexpr("_",row.names(all_psi), fixed=TRUE)))
      lambda <- plyr::ddply(all_psi, "unit", function (x){
          # sort by time
          x <- x[ order(x$time),]
          # estimate the ratio 
          est <- x$est[-1]/x$est[-nrow(x)]
          # we cannot estimate the standard error because we don't know the covariance between the estimates of psi for a unit
          lambda <- data.frame(
                         time = x$time[-nrow(x)],
                         est= est,
                         se = NA,
                         lower_0.95=NA,
                         upper_0.95=NA,
                         stringsAsFactors=FALSE)
          lambda
      })
      row.names(lambda) <- paste(lambda$unit,"_",lambda$time,sep="")
      model$derived$lambda <- lambda
      
      model$derived$lambda$unit <- NULL
      model$derived$lambda$time <- NULL
      
      # Estimate the chance in the occupancy rate on the logit scale (lamba') for each site
      logit <- function(x) log(x/(1-x))
      expit <- function(x) 1/(1+exp(-x))
      
      lambdap <- plyr::ddply(all_psi, "unit", function (x){
        # sort by time
        x <- x[ order(x$time),]
        # estimate the ratio 
        est <-exp(logit(x$est[-1])-logit(x$est[-nrow(x)]))
        # we cannot estimate the standard error because we don't know the covariance between the estimates of psi for a unit
        lambdap <- data.frame(
          time = x$time[-nrow(x)],
          est= est,
          se = NA,
          lower_0.95=NA,
          upper_0.95=NA,
          stringsAsFactors=FALSE)
        lambdap
      })
      #browser()
      row.names(lambdap) <- paste(lambdap$unit,"_",lambdap$time,sep="")
      model$derived$lambdap <- lambdap
      model$derived$lambdap$unit <- NULL
      model$derived$lambdap$time <- NULL
  }
  
  if("do4" %in% class(model)){ # single species multi season random occupancy model
    # The real psi has all of the estimates already.
    model$derived$all_psi <- model$real$psi
    
    # Estimate the occupancy growth rate (lambda) for each site and year
    all_psi <- model$derived$all_psi
    all_psi$unit <- substr(row.names(all_psi),1,-1+regexpr("_",row.names(all_psi), fixed=TRUE))
    all_psi$time <- as.numeric(substring(row.names(all_psi),1+regexpr("_",row.names(all_psi), fixed=TRUE)))
    lambda <- plyr::ddply(all_psi, "unit", function (x){
      # sort by time
      x <- x[ order(x$time),]
      # estimate the ratio 
      est <- x$est[-1]/x$est[-nrow(x)]
      # we cannot estimate the standard error because we don't know the covariance between the estimates of psi for a unit
      lambda <- data.frame(
        time = x$time[-nrow(x)],
        est= est,
        se = NA,
        lower_0.95=NA,
        upper_0.95=NA,
        stringsAsFactors=FALSE)
      lambda
    })
    row.names(lambda) <- paste(lambda$unit,"_",lambda$time,sep="")
    model$derived$lambda <- lambda
    
    model$derived$lambda$unit <- NULL
    model$derived$lambda$time <- NULL
    
    # Estimate the chance in the occupancy rate on the logit scale (lamba') for each site
    logit <- function(x) log(x/(1-x))
    expit <- function(x) 1/(1+exp(-x))
    
    lambdap <- plyr::ddply(all_psi, "unit", function (x){
      # sort by time
      x <- x[ order(x$time),]
      # estimate the ratio 
      est <-exp(logit(x$est[-1])-logit(x$est[-nrow(x)]))
      # we cannot estimate the standard error because we don't know the covariance between the estimates of psi for a unit
      lambdap <- data.frame(
        time = x$time[-nrow(x)],
        est= est,
        se = NA,
        lower_0.95=NA,
        upper_0.95=NA,
        stringsAsFactors=FALSE)
      lambdap
    })
    #browser()
    row.names(lambdap) <- paste(lambdap$unit,"_",lambdap$time,sep="")
    model$derived$lambdap <- lambdap
    model$derived$lambdap$unit <- NULL
    model$derived$lambdap$time <- NULL
  }
  
  model
}


RPresence.modAvg.derived <- function(aic.set, param, estimate="est", se="se"){
   # compute model average of derived parameters
   require(plyr)
   # Extract the estimates and standard errors
   est     <- plyr::laply(aic.set$models, function(x){ x$derived[[param]][,estimate]})
   est.se  <- plyr::laply(aic.set$models, function(x){ x$derived[[param]][,se      ]})
   wgt <- plyr::laply(aic.set$models, function(x,aic.set){
       aic.set$table$wgt[ which(x$modname == aic.set$table$Model) ]
   }, aic.set=aic.set)
   rnames  <- row.names(aic.set$models[[1]]$derived[[param]])
   
   # Compute the model average
   ma_est <-wgt %*% as.matrix(est)
   var2   <- (est - as.matrix(ma_est)[rep(1,length(wgt)),])^2
   ma_se  <-wgt %*% sqrt(est.se^2 + var2)
   res <- data.frame(est=as.vector(ma_est),
                     se =as.vector(ma_se),
                     lower_0.95=as.vector(ma_est - qnorm(0.975)*ma_se), 
                     upper_0.95=as.vector(ma_est + qnorm(0.975)*ma_se),
                     stringsAsFactors=FALSE)
   row.names(res) <- rnames
   res
    }
   
