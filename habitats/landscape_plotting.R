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

setwd("~/Documents/GitHub/thesis_simulation")
set.seed(2024)

low.hab <- readRDS("habitats/low_hab.rds")
med.hab <- readRDS("habitats/med_hab.rds")
high.hab <- readRDS("habitats/high_hab.rds")

n=25

#~~~~~~~~
# Examine the objects
#~~~~~~~~
# par(mfrow=c(1, 3))
# plot(low.hab)
# 
# plot(medium.hab)
# 
# plot(high.hab)
# par(mfrow=c(1,1))
# 
# summary(values(low.hab))
# summary(values(medium.hab))
# summary(values(high.hab))
# 
# par(mfrow=c(1,3))
# plot(low.hab, col=rev(rainbow(99, start=0,end=1)),
#       breaks=seq(min(values(high.hab)),max(values(high.hab)), length.out=100))
# plot(medium.hab, col=rev(rainbow(99, start=0,end=1)),
#      breaks=seq(min(values(high.hab)),max(values(high.hab)), length.out=100))
# plot(high.hab, col=rev(rainbow(99, start=0,end=1)),
#      breaks=seq(min(values(high.hab)),max(values(high.hab)), length.out=100))
# par(mfrow=c(1,1))
vals.sim <- c(values(low.hab), values(med.hab), values(high.hab))
# cols <- sequential_hcl(length(unique(vals.sim)))



low.hab2 <- low.hab
# low.hab2$cols <- vals.sim[1:2500]
low.hab2$cols <- values(low.hab)

l <- ggplot() +
  geom_spatraster(data=low.hab2, aes(fill=cols)) +
  scale_fill_continuous_divergingx(palette="Geyser",
                               breaks=round(seq(15, -15, length.out=7)),
                               limits=c(min(vals.sim), max(vals.sim)), rev=T) +
  theme_classic() +
  labs(title="Low") +
  lims(x=c(-n, n), y=c(-n, n)) +
  guides(fill=guide_legend(title="Value")) +
  theme(legend.position = "none")

png("plots/low_var_hab.png", 6, 6, units="in", res=600)
print(l)
dev.off()


medium.hab2 <- med.hab
medium.hab2$cols <- values(med.hab)
m <- ggplot() +
  geom_spatraster(data=medium.hab2, aes(fill=cols)) +
  scale_fill_continuous_divergingx(palette="Geyser",
                               breaks=round(seq(15, -15, length.out=7)),
                               limits=c(min(vals.sim), max(vals.sim)), rev=T) +
  theme_classic() +
  labs(title="Medium") +
  lims(x=c(-n, n), y=c(-n, n)) +
  guides(fill=guide_legend(title="Value")) +
  theme(legend.position = "none")

png("plots/med_var_hab.png", 6, 6, units="in", res=600)
print(m)
dev.off()


high.hab2 <- high.hab
high.hab2$cols <- vals.sim[5001:7500]
h <- ggplot() +
  geom_spatraster(data=high.hab2, aes(fill=cols)) +
  scale_fill_continuous_divergingx(palette="Geyser",
                               breaks=round(seq(15, -15, length.out=7)),
                               limits=c(min(vals.sim), max(vals.sim)), rev=T) +
  theme_classic() +
  # labs(title="High variation habitat") +
  labs(title="High") +
  lims(x=c(-n, n), y=c(-n, n)) +
  guides(fill=guide_legend(title="Value")) +
  theme(legend.position = "bottom")
 # theme(legend.position = "none")

png("plots/legend.png", 6, 6, units="in", res=600)
print(h)
dev.off()
