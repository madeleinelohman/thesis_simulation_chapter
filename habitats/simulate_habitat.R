#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Libraries and set up
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(ClustGeo)
library(clusterCrit)
library(colorspace)
library(cowplot)
library(dbscan)
library(geoR)
library(MASS)
library(mvtnorm)
library(raster)
library(rgeoda)
library(sf)
library(spdep)
library(tidyverse)
library(tidyterra)
library(terra)

setwd("/Users/madelienelohman/Documents/GitHub/thesis_simulation")
set.seed(2048)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulate habitat
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

n = 25 # Grid dimensions are n*2 x n*2
phi = 5 # Range
sig.levels = seq(5, 25, length.out=3) # Partial sill/variance
sig.names = c("low", "medium", "high") # Names for the level of variance

for(i in 1:length(sig.levels)){
  ### Create an empty raster
  r=rast(ncol=n*2,nrow=n*2,extent=c(-n,n,-n,n))
  r[]=0
  s.r=crds(r)
  ### Fill it using a Gaussian random field
  rf=grf(n, xyFromCell(r,1:ncell(r)), cov.pars=c(sig.levels[i], phi), nugget=3,
         cov.model="gaussian")
  ### Populate a new matrix and create an object out of it for each level of variation
  hab=r
  hab[]=rf$data
  assign(paste0(sig.names[i],".hab"), hab)
}

saveRDS(low.hab, "low_hab.rds")
saveRDS(medium.hab, "med_hab.rds")
saveRDS(high.hab, "high_hab.rds")

