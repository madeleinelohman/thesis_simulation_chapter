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

load("scoring_res.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set up scoring dataframe
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~
# Clean data
#~~~~~~~~~~~
### Create appropriate factor levels so things show up in plots correctly
ints$var.level <- factor(ints$var.level, levels=c("high", "medium", "low"))
ints$alg.level <- paste(ints$Algorithm, ints$var.level)
ints$alg.level <- factor(ints$alg.level, levels=rev(c("Grid low", "Grid medium", "Grid high",
                                                  "Voronoi low", "Voronoi medium", "Voronoi high",
                                                  "AZP low", "AZP medium", "AZP high",
                                                  "SKATER low", "SKATER medium", "SKATER high",
                                                  "REDCAP low", "REDCAP medium", "REDCAP high")))

### Clean up names of optimization indicies
ints$Index[which(ints$Index == "Calinski_Harabasz")] <- "Calinski-Harabasz"
ints$Index[which(ints$Index == "Davies_Bouldin")] <- "Davies-Bouldin"
ints$Index[which(ints$Index == "SD_Scat")] <- "SD"


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plotting
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~
# Total number of clusters
#~~~~~~~~~~
tot.clust <- ints %>%
  ggplot(aes(x=Index, y=alg.level, fill=Clusters)) +
  geom_tile() +
  scale_fill_continuous_divergingx(palette="Geyser", mid=25) +
  labs(y="Algorithm/Variation level combination", title="Number of clusters") +
  theme_classic()
  
png("plots/all_alg_clusts_ng.png", 7.5, 6, units="in", res=600)
print(tot.clust)
dev.off()

#~~~~~~~~~~
# Goodness of fit
#~~~~~~~~~~
rs.heat <- ggplot(ints, aes(x=Index, y=alg.level, fill=RS)) + 
  geom_tile() +
  scale_fill_continuous_sequential(palette="Mako") +
  labs(y="Algorithm/Variation level combination", title="Goodness of fit (RS)") +
  theme_classic()

png("plots/rs_heatmap.png", 7.5, 7, units="in", res=600)
print(rs.heat)
dev.off()


#~~~~~~~~~~
# Means for different categories
#~~~~~~~~~~
### Mean RS score by index
mus.scores.rs <- ints %>%
  group_by(Index) %>%
  summarise(mu_RS=round(mean(RS), 3))

### Meancluster number by index (excluding grid)
mus.scores.clust <- ints %>%
  group_by(Index) %>%
  summarise(mu_clust=mean(Clusters))

### Mean RS scores and cluster numbers by variation level (excluding grid)
levels.scores <- ints %>%
  group_by(var.level) %>%
  summarise(mu_clust=mean(Clusters), mu_RS=mean(RS))



#~~~~~~~~~~
# Individual algorithm details
#~~~~~~~~~~
ints %>%
  filter(Algorithm == 'SKATER') %>%
  group_by(var.level) %>%
  summarize(mu_clust=mean(Clusters), sd_clust=sd(Clusters))

ints %>%
  filter(Algorithm == 'REDCAP') %>%
  group_by(var.level) %>%
  summarize(mu_clust=mean(Clusters), sd_clust=sd(Clusters))


ints %>%
  filter(Algorithm == 'AZP') %>%
  group_by(var.level) %>%
  summarize(mu_clust=mean(Clusters), sd_clust=sd(Clusters))

ints %>%
  filter(Algorithm == 'Voronoi') %>%
  group_by(var.level) %>%
  summarize(mu_clust=mean(Clusters), sd_clust=sd(Clusters))
ints %>%
  filter(Algorithm == 'Voronoi') %>%
  summarize(mean(RS))

ints %>%
  filter(Algorithm == 'Grid') %>%
  group_by(var.level) %>%
  summarize(mu_clust=mean(Clusters), sd_clust=sd(Clusters))


#~~~~~~~~~~
# Individual index details
#~~~~~~~~~~
ch <- ints %>%
  filter(Index=="Dunn") %>%
  group_by(Algorithm)


