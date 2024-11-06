

v.n.clust <- function(x, min.clust, max.clust, index, n.iter){
  
  v <- as.polygons(x, aggregate=F) 
  s <- sf::st_as_sf(v)
  x1 <- matrix(values(x), nrow=50, ncol=50, byrow=T)
  
  ### Find point values from the centroid of each cell
  cents.pts <- crds(x)
  cents <- terra::extract(x, cents.pts)
  cents <- cbind(cents, cents.pts)
  cents <- st_as_sf(cents, coords=c('x', 'y'))
  
  ### Get input information
  queen_w <- queen_weights(s) # Queen weighted neighborhood
  data <- s[c('lyr.1')] # Data of interest
  
  env <- list(matrix(c(-25,-25,25,-25,25,25,-25,25, -25,-25), ncol=2, byrow=T))
  bbox = st_sfc(st_polygon(env))
  
  eval <- data.frame(clusters=min.clust:max.clust, Index=NA)
  
  for(i in min.clust:max.clust){
    score <- rep(NA, times=n.iter)
    for(k in 1:n.iter){
      samps <- terra::spatSample(x, size=i, method="random", as.df=T, xy=T)
      samps <- st_as_sf(samps, coords=c('x', 'y'))
      
      v <- st_combine(st_geometry(samps)) %>%
        st_voronoi(envelope=bbox) %>%
        st_collection_extract()
      #v <- v[unlist(st_intersects(samps, v))]
      v <- st_as_sf(v)
      
      matchy <- st_contains(v, cents)
      touchy <- st_touches(v, cents)
      clust <- matrix(NA, n*2, n*2)
      for(j in 1:i){
        want.these <- c(matchy[[j]], touchy[[j]])
        clust[want.these] <- j
      }
      
      new.clusts <- order.clusts(c(clust))
      if(index == "Dunn"){
        score[k] <- dunn.score(x, new.clusts$new)
      }else{
        score[k] <- unlist(intCriteria(x1, new.clusts$new, crit=index))
      }
    }
    if(all(is.nan(score))){
      eval$Index[i-(min.clust-1)] <- NA
    }else{
      score <- score[which(!is.nan(score))]
      eval$Index[i-(min.clust-1)] <- mean(score, na.rm=T)
    }
  }
  
  if(index == "Dunn"){
    best.want <- which.max(eval$Index == max(eval$Index))
  }else{
    best.want <- which.max(eval$Index == eval$Index[bestCriterion(eval$Index, index)])
  }
  
  n.clust.want <- eval[best.want, "clusters"]
  
  return(n.clust.want)
}
