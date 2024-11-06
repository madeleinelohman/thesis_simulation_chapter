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
set.seed(2024)


source("algs/grd_alg.R")
source("algs/voronoi_alg.R")
source("algs/azp_alg.R")
source("algs/SKATER_alg.R")
source("algs/REDCAP_alg.R")
source("misc_funcs/order_clusts.R")
source("misc_funcs/dunn_score.R")

### Read in simulated rasters
low.hab <- readRDS("habitats/low_hab.rds")
medium.hab <- readRDS("habitats/med_hab.rds")
high.hab <- readRDS("habitats/high_hab.rds")


n = dim(low.hab)[1]/2
min.clust = 3
max.clust = 50

start.time <- Sys.time()
indices <- c("Calinski_Harabasz", "Davies_Bouldin", "Dunn", "SD_Scat", "Silhouette")

for(k in 1:length(indices)){
  
  index=indices[k]


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Grid
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  low.grid <- grid.alg(low.hab, "Low", n, index, min.clust, max.clust)
  low.grid$score
  
  medium.grid <- grid.alg(medium.hab, "Medium", n, index, min.clust, max.clust)
  medium.grid$score
  
  high.grid <- grid.alg(high.hab, "High", n, index, min.clust, max.clust)
  high.grid$score
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Voronoi tesselation
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  n.iter=25
  
  low.v <- v.alg(low.hab, "Low", min.clust, max.clust, n, index, n.iter)
  low.v$score
  
  medium.v <- v.alg(medium.hab, "Medium", min.clust, max.clust, n, index, n.iter)
  medium.v$score
  
  high.v <- v.alg(high.hab, "High", min.clust, max.clust, n, index, n.iter)
  high.v$score
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # AZP
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  low.azp <- azp.alg(low.hab, "Low", min.clust, max.clust, index)
  low.azp$score
  
  medium.azp <- azp.alg(medium.hab, "Medium", min.clust, max.clust, index)
  medium.azp$score
  
  high.azp <- azp.alg(high.hab, "High", min.clust, max.clust, index)
  high.azp$score
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # SKATER
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  low.skater <- skater.alg(low.hab, "Low", min.clust, max.clust, index, n)
  low.skater$score
  low.skater$n.clust
  
  medium.skater <- skater.alg(medium.hab, "Medium", min.clust, max.clust, index, n)
  medium.skater$score
  medium.skater$n.clust
  
  high.skater <- skater.alg(high.hab, "High", min.clust, max.clust, index, n)
  high.skater$score
  high.skater$n.clust
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # REDCAP
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  low.red <- redcap.alg(low.hab, "Low", min.clust, max.clust, index, n)
  low.red$score
  
  medium.red <- redcap.alg(medium.hab, "Medium", min.clust, max.clust, index, n)
  medium.red$score
  
  high.red <- redcap.alg(high.hab, "High", min.clust, max.clust, index, n)
  high.red$score
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Save everything
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  save.image(paste0("cluster_results/cluster_res_",index,".RData"))
  
  
  end.time <- Sys.time()
  end.time - start.time

}