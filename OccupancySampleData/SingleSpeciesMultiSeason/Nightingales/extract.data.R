# extract the data on nightinggale from the blmeco package

library(blmeco)
library(plyr)

data("nightingales")

# This is a 3 dimensional array {sites, year, visits (1-8)}

night.hist <- adply(nightingales, 1, function(x){
      dhist <- as.vector(unlist(t(x)))
      dhist
})
head(night.hist)

night.hist <- plyr::rename(night.hist, c("X1"="Site"))
head(night.hist)
nightingales[1,,]


write.csv(night.hist, 'nightingale.csv', row.names=FALSE)