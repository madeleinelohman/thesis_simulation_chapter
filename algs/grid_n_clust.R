

grid.n.clust <- function(x, index, n, min.clust, max.clust){
  ### Set up values for for loop
  rast.vals <- matrix(values(x), n*2, n*2) ### Make a matrix out of the raster values
  total.mean <- mean(values(x)) ### Find the mean of the raster
  
  ### What are potential window sizes?
  clusts <- min.clust:max.clust
  # clusts <- clusts[clusts %% 2 == 0 | clusts %% 5 == 0]
  clusts <- clusts[clusts %% 2 == 0]
  
  eval <- data.frame(clusters=clusts, Index=NA)
  
  for(k in 1:length(clusts)){
    ### Set up matrices for for loop
    grp.lvl <- round(sqrt((n*2)^2/clusts[k]))
    splts <- split(1:(n*2), ceiling(seq_along(1:(n*2))/grp.lvl))
    n.grps <- length(splts)
    
    lats <- longs <- unlist(lapply(splts, function(x) first(x)))
    lats[n.grps+1] <- longs[n.grps+1] <- 51
    
    ### Set up matrices for the for loop
    grd <- matrix(NA, n.grps, n.grps) # New grid
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
      }
    }
    
    new.clusts <- order.clusts(c(clusters))
    
    if(index == "Dunn"){
      eval$Index[k-(min.clust-1)] <- dunn.score(x, new.clusts$new)
    }else{
      eval$Index[k-(min.clust-1)] <- unlist(intCriteria(rast.vals, new.clusts$new, crit=index)) 
    }
  }
  
  if(index == "Dunn"){
    best.want <- which.max(eval$Index == max(eval$Index, na.rm=T) &
                             !is.nan(eval$Index))
  }else{
    best.want <- which.max(eval$Index == eval$Index[bestCriterion(eval$Index, index)] &
                             !is.nan(eval$Index))
  }
  
  n.clust.want <- eval[best.want, "clusters"]
  
  return(n.clust.want)
}