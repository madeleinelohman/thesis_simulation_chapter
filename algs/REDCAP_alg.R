redcap.alg <- function(x, level, min.clust, max.clust, index, n){
  source("algs/redcap_n_clust.R")
  #~~~~~~~~~~~~~~~~~~~~
  # Decide number of clusters
  #~~~~~~~~~~~~~~~~~~~~
  n.clust <- redcap.n.clust(x, min.clust, max.clust, index)
  
  #~~~~~~~~~~~~~~~~~~~~
  # Run method
  #~~~~~~~~~~~~~~~~~~~~
  ### Transform raster to sf class
  # Need to go through sp class to get to sf class from SpatRaster 
  v <- as.polygons(x, aggregate=F) 
  s <- sf::st_as_sf(v)
  
  ### Get input information
  queen_w <- queen_weights(s) # Queen weighted neighborhood
  data <- s[c('lyr.1')] # Data of interest
  
  ### Run REDCAP!
  cr <- redcap(n.clust, queen_w, data, "fullorder-completelinkage") 
  
  #~~~~~~~~~~~~~~~~~~~~
  # Plotting and plotting setup
  #~~~~~~~~~~~~~~~~~~~~
  ### Create a dataframe with old values of raster and new clusters
  pred.df <- data.frame(org=values(x), new=cr$Clusters)
  names(pred.df)[1] <- "org"
  
  ### Summarize these clusters by their mean values
  summ.pred <- pred.df %>%
    group_by(new) %>%
    summarize(vals = mean(org))
  
  ### Associate this new mean with the clusters
  pred.df$new.mean <- NA
  clusts <- unique(summ.pred$new)
  for(i in 1:length(clusts)){
    pred.df[which(pred.df$new == clusts[i]), 3] <- summ.pred[which(summ.pred$new == clusts[i]), 2]
  }
  
  ### Create a new raster with these new means for each cluster
  new.hab <- x
  values(new.hab) <- pred.df$new.mean
  
  ### Plot!
  par(mfrow=c(1,2))
  plot(x, main=paste("Simulated habitat:", level))
  plot(new.hab, main=paste("Predicted clusters:", level))
  par(mfrow=c(1,1))
  
  # png("plots/redcap_pred.png", 6, 6, "in", res=600)
  # plot(hab.cat4, main="REDCAP predicted habitat")
  # dev.off()
  
  p <- ggplot() +
    geom_spatraster(data=new.hab) +
    scale_fill_binned_divergingx(palette="Geyser") +
    theme_classic() + 
    labs(title=paste(level, "variation habitat"), subtitle="REDCAP") +
    lims(x=c(-n, n), y=c(-n, n)) +
    guides(fill=guide_legend(title="Value"))
  
  png(paste0("plots/redcap_plots/redcap_",level,"_",index,".png"), 6, 6, res=600, units="in")
  print(p)
  dev.off()
  
  #~~~~~~~~~~~~~~~~~~~~
  # Validate method
  #~~~~~~~~~~~~~~~~~~~~
  ### The following code is similar to the grid algorithm validation
  rast.vals <- matrix(values(x), n*2, n*2)
  total.mean <- mean(rast.vals)
  pred <- matrix(cr$Clusters, nrow=n*2)
  
  ssw <- NA
  ssb <- NA
  for(i in 1:n.clust){
    ind.want <- which(pred == i, arr.ind=T)
    m.want <- mean(rast.vals[ind.want])
    
    ssw[i] <- sum((rast.vals[ind.want] - m.want)^2)
    ssb[i] <- (m.want - total.mean)^2 * nrow(ind.want)
  }
  ss.total <- sum(ssb + ssw)
  ssb.sst <- sum(ssb) / ss.total
  
  new.clusts <- order.clusts(cr$Clusters)
  if(index == "Dunn"){
    score <- dunn.score(x, new.clusts$new)
  }else{
    score <- unlist(intCriteria(x1, new.clusts$new, crit=index)) 
  }
  
  ### Return objects
  return(list(new.hab=pred, ssw=ssw, ssb=ssb, sst=ss.total, ssb.sst=ssb.sst,
              clusters=cr$Clusters, n.clust=n.clust, p=p, score=score))
}
