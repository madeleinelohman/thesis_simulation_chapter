grid.alg <- function(x, level, n, index, min.clust, max.clust){
  
  #~~~~~~~~~~~~~~~~~~~~
  # Decide number of clusters
  #~~~~~~~~~~~~~~~~~~~~
  source("algs/grid_n_clust.R")
  n.clust <- grid.n.clust(x, index, n, min.clust, max.clust)
  
  #~~~~~~~~~~~~~~~~~~~~
  # Run method
  #~~~~~~~~~~~~~~~~~~~~
  rast.vals <- matrix(values(x), n*2, n*2) ### Make a matrix out of the raster values
  total.mean <- mean(values(x)) ### Find the mean of the raster

  grp.lvl <- round(sqrt((n*2)^2/n.clust))
  splts <- split(1:(n*2), ceiling(seq_along(1:(n*2))/grp.lvl))
  n.grps <- length(splts)
  
  lats <- longs <- unlist(lapply(splts, function(x) first(x)))
  lats[n.grps+1] <- longs[n.grps+1] <- 51
  
  ### Set up matrices for the for loop
  # grp.lvl = Window to group cells by (new grid dimension sizes)
  grd <- matrix(NA, n.grps, n.grps) # New grid
  ssw <- matrix(NA, n.grps, n.grps) # Within-cluster sum of squares
  ssb <- matrix(NA, n.grps, n.grps) # Between-cluster sum of squares
  clusters <- matrix(0, n*2, n*2) # Cluster assignment
  
  
  for(i in 1:(length(lats)-1)){
    for(j in 1:(length(lats)-1)){
      ### Identify the cell ranges to aggregate
      ind.lat <- lats[i]:(lats[i+1]-1) 
      ind.long <- longs[j]:(longs[j+1]-1)
      
      ### Create grid
      mean.local.vals <- mean(rast.vals[ind.lat, ind.long]) # Find the mean from the original raster
      grd[i,j] <- mean.local.vals # Place the mean in the new grid
      clusters[ind.lat, ind.long] <- max(clusters) + 1 # Assign a cluster number
      
      ### Validate method
      ssw[i,j] <- sum((rast.vals[ind.lat, ind.long] - mean.local.vals)^2) # Calculate individual cell SSW
      ssb[i,j] <- (mean.local.vals - total.mean)^2 * grp.lvl^2 # Calculate individual cell SSB
    }
  }
  ss.total <- sum(ssb + ssw) # Total within-cluster sum of squares
  ssb.sst <- sum(ssb) / ss.total # The ratio of between to total sum of squares
  
  #~~~~~~~~~~~~~~~~~~~~
  # Rasterize results
  #~~~~~~~~~~~~~~~~~~~~
  new.hab <- rast(grd,extent=c(-n,n,-n,n))
  
  p <- ggplot() +
    geom_spatraster(data=new.hab) +
    scale_fill_binned_divergingx(palette="Geyser") +
    theme_classic() + 
    labs(title=paste(level, "variation habitat"), subtitle="Grid") +
    lims(x=c(-n, n), y=c(-(n+1), n)) +
    guides(fill=guide_legend(title="Value"))
  
  png(paste0("plots/grid_plots/grid_",level,index,".png"), 6, 6, res=600, units="in")
  print(p)
  dev.off()
  
  
  new.clusts <- order.clusts(c(clusters))
  if(index == "Dunn"){
    score <- dunn.score(x, new.clusts$new)
  }else{
    score <- unlist(intCriteria(rast.vals, new.clusts$new, crit=index)) 
  }
  
  ### Return objects
  return(list(new.hab=new.hab, ssw=ssw, ssb=ssb, sst=ss.total, ssb.sst=ssb.sst,
              clusters=clusters, p=p, score=score, n.clust=max(clusters)))
}
