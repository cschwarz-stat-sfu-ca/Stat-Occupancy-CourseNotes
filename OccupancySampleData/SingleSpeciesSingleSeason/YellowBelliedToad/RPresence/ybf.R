# Single Species Single Season Occupancy 

# Yellow-bellied toad
# Single Season Single Season occupancy 

#  RPresence package

# 2018-11-15 Code contributed by Carl James Schwarz (cschwarz.stat.sfu.cs@gmail.com)

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------



library(readxl)
library(RPresence)
library(ggplot2)

# get the data read in
# Data for detections should be a data frame with rows corresponding to sites
# and columns to visits.
# The usual 1=detected; 0=not detected; NA=not visited is used.

input.history <- readxl::read_excel(file.path("..","YellowBelliedToad.xlsx"),
                                    sheet="detections")

# do some basic checks on your data 
# e.g. check number of sites; number of visits etc
nrow(input.history)
ncol(input.history)
range(input.history, na.rm=TRUE)
sum(is.na(input.history))

head(input.history)

# Get the site level covariates - none.



# Get the site x visit (sampling) covariates. 
jdate <- readxl::read_excel(file.path("..","YellowBelliedToad.xlsx"),
                                    sheet="JulianDate")
range(jdate)

# Observational covariate needs to be "stacked" so that sites1...siteS for survey occastion 1
# are then followed by covariate at survey occastion 2 for sites1...siteS, etc

survey.cov <- data.frame(site =rep(1:nrow(input.history) , ncol(input.history)),
                         visit =as.character(rep(1:ncol(input.history), each=nrow(input.history))),  # notice we make a character 
                         jdate =as.vector(unlist(jdate)), 
                         stringsAsFactors=FALSE)
head(survey.cov[,c("visit","site","jdate")],n=50)
str(survey.cov)

# standarize jdate
survey.cov$jdateS <- (survey.cov$jdate - 150)/10

# check that missing values in history and observer covariates align
select <- is.na(as.vector(unlist(input.history)))
survey.cov[select,]
sum(is.na(survey.cov[!select,]))


# Create the *.pao file
ybf.pao <- RPresence::createPao(input.history,
                                 #unitcov=site_covar, # no site covariates
                                 survcov=survey.cov, 
                                 title='ybf SSSS')
ybf.pao

#--------------------------------------------------------------------------
# Fit some models.
# Note that formula DO NOT HAVE AN = SIGN
mod.pdot <- RPresence::occMod(model=list(psi~1, p~1), type="so", data=ybf.pao)
summary(mod.pdot)

# look at estimated occupancy probability. RPresence gives for EACH site in case it depends on covariates
mod.pdot.psi <-mod.pdot$real$psi[1,]  # occupancy probability
mod.pdot.psi


# look at the estimated probability of detection. It gives an estimate for every site at very visit
mod.pdot.p   <- mod.pdot$real$p[seq(1, by=nrow(input.history), length.out=ncol(input.history)),]
mod.pdot.p

# alternatively
RPresence::print_one_site_estimates(mod.pdot, site = 1)


#-------
# Model where p(t) varies across survey occasions  
# 
mod.pt.psidot <- RPresence::occMod(model=list(psi~1, p~SURVEY), type="so", data=ybf.pao)
summary(mod.pt.psidot)

mod.pt.psidot$real$psi[1:5,]
mod.pt.psidot$real$p[seq(1, by=nrow(input.history), length.out=ncol(input.history)),]

print_one_site_estimates(mod.pt.psidot, site = 1)
fitted(mod.pt.psidot, param="psi")[1,]

#-------
# Model where p varies by julian date in a linear form 
# 
mod.pjdate.psidot <- RPresence::occMod(model=list(psi~1, p~jdateS), type="so", data=ybf.pao)
summary(mod.pjdate.psidot)

mod.pjdate.psidot$real$psi[1:5,]
mod.pjdate.psidot$real$p[seq(1, by=nrow(input.history), length.out=ncol(input.history)),]

print_one_site_estimates(mod.pjdate.psidot, site = 1)
fitted(mod.pjdate.psidot, param="psi")[1:5,]



#-------
# Model where p varies by julian date in a quadrat form 
# Notice we need the I() around the square term in the model to force R to evaluate it properly.
# See https://stackoverflow.com/questions/8055508/in-r-formulas-why-do-i-have-to-use-the-i-function-on-power-terms-like-y-i
mod.pjdate2.psidot <- RPresence::occMod(model=list(
                                         psi~1,
                                         p~jdateS + I(jdateS^2)), 
                                        type="so", data=ybf.pao)
summary(mod.pjdate2.psidot)

mod.pjdate2.psidot$beta

mod.pjdate2.psidot$real$psi[1:5,]
mod.pjdate2.psidot$real$p[seq(1, by=nrow(input.history), length.out=ncol(input.history)),]

print_one_site_estimates(mod.pjdate2.psidot, site = 1)
fitted(mod.pjdate2.psidot, param="psi")[1:5,]





#------
# Model averaging
models<-list(mod.pdot,
             mod.pjdate.psidot, 
             mod.pjdate2.psidot, 
             mod.pt.psidot)

results<-RPresence::createAicTable(models)
summary(results)

RPresence::modAvg(results, param="psi")[1:5,]

# look at detectability 
ma.p <- RPresence::modAvg(results, param="p")

# we need to extract the site and survey and merge with survey.cov to get jdate
head(ma.p)
ma.p$site <- as.numeric( substring(row.names(ma.p), 4+regexpr("unit", row.names(ma.p))))
ma.p$visit<- substr(row.names(ma.p),2,2)
head(ma.p)

ma.p <- merge(ma.p, survey.cov)
head(ma.p)

ggplot(data=ma.p, aes(x=jdate, y=est))+
  ggtitle("Estimated detection as a function of julian date")+
  geom_point()+
  geom_ribbon(aes(ymin=lower_0.95, ymax=upper_0.95), alpha=0.2)+
  ylim(0,1)

