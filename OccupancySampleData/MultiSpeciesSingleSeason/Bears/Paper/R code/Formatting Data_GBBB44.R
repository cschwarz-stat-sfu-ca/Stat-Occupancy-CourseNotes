setwd("C:/")

# importing detection / non-detection data
grizzly <- as.matrix(read.csv('four_day_grizzly.csv', F))
black <- as.matrix(read.csv('four_day_black.csv', F))
motorised <- as.matrix(read.csv('four_day_motorized.csv', F))
non_motorised <- as.matrix(read.csv('four_day_non_motorized.csv', F))

# total number of camera days
cday <- apply(grizzly, 1, function(x) sum(is.na(x) == F))

# importing covariate data
psi.cov <- read.csv('covariates.csv')
# p.cov <- read.csv('p covariates.csv')

# CREATING DESIGN MATRIX
X.array <- array(0, dim = c(nrow(grizzly), 16, 15))

droad <- log(psi.cov$droad+1)
dstream <- log(psi.cov$dstream+1)
elev <- scale(psi.cov$elev)[, 1]
tpi <- scale(psi.cov$tpi)[, 1]
NDVI <- scale(psi.cov$NDVI)[, 1]
prot<-psi.cov$Land_type

# importing the general form of the design matrix
library(xlsx)
dm <- read.xlsx('Design Matrix_GBBB.xlsx', 2, rowIndex = 3:18, colIndex = 2:16,
                header = F)

# filling in the design matrix
for(i in 1:nrow(grizzly)){
  for(j in 1:15){
    X.array[i, j, dm[j, ] == 1] <-
      c(1, droad[i], dstream[i], elev[i], NDVI[i],
        1, droad[i], dstream[i], elev[i], NDVI[i],
        1, prot[i],
        1, prot[i],
        1
        
      )[dm[j, ] == 1]
  }
}

# OBSERVED Y
y1max <- apply(grizzly, 1, max, na.rm = T)
y2max <- apply(black, 1, max, na.rm = T)
y3max <- apply(motorised, 1, max, na.rm = T)
y4max <- apply(non_motorised, 1, max, na.rm = T)

# Vectorizing detection / non-detection
Y1 <- Y2 <- Y3 <- Y4 <- prot <- numeric()
for(i in 1:nrow(grizzly)){
  Y1 <- c(Y1, grizzly[i, 1:cday[i]])
  Y2 <- c(Y2, black[i, 1:cday[i]])
  Y3 <- c(Y3, motorised[i, 1:cday[i]])
  Y4 <- c(Y4, non_motorised[i, 1:cday[i]])
  prot = c(prot, rep(psi.cov$Land_type[i], cday[i]))
}

# Stan does not support ragged indexing.  This is a work-around
start <- 1
for(i in 2:nrow(grizzly)){
  start <- c(start, sum(cday[1:(i - 1)]) + 1)
}

data <- list(
  K = ncol(X.array[1, , ]),
  L = 2,
  N = nrow(grizzly),
  NJ = length(Y1),
  S = 16,
  obs = cday,
  start = start,
  x = X.array,
  I1 = y1max, I2 = y2max, I3 = y3max, I4 = y4max,
  Y1 = Y1, Y2 = Y2, Y3 = Y3, Y4 = Y4,
  prot=prot
  )

inits <- function(){
  list(
    a1 = array(rnorm(1), 1),
    a2 = array(rnorm(1), 1),
    a3 = rnorm(2),
    a4 = rnorm(2),
    a1_c1 = rnorm(1),
    a1_c2 = rnorm(1),
    a2_c1 = rnorm(1),
    a2_c2 = rnorm(1),
    a2_c3 = rnorm(1),
    beta = rnorm(ncol(X.array[1, , ])))
}

params <- c('a1', 'a2', 'a3', 'a4', 'a1_c1', 'a1_c2', 'a2_c1', 'a2_c2', 'a2_c3', 'beta', 'll')
