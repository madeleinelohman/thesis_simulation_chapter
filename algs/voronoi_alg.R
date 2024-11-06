v.alg <- function(x, level, min.clust, max.clust, n, index, n.iter){
  #~~~~~~~~~~~~~~~~~~~~
  # Decide number of points
  #~~~~~~~~~~~~~~~~~~~~
  source("algs/voronoi_n_clust.R")
  n.clust <- v.n.clust(x, min.clust, max.clust, index, n.iter)
  v <- as.polygons(x, aggregate=F) 
  s <- sf::st_as_sf(v)
  
  #~~~~~~~~~~~~~~~~~~~~
  # Run method
  #~~~~~~~~~~~~~~~~~~~~
  samps <- terra::spatSample(x, size=n.clust, method="random", as.df=T, xy=T)
  samps <- st_as_sf(samps, coords=c('x', 'y'))
  
  env <- list(matrix(c(-25,-25,25,-25,25,25,-25,25, -25,-25), ncol=2, byrow=T))
  bbox <- st_sfc(st_polygon(env))
  
  v <- st_combine(st_geometry(samps)) %>%
    st_voronoi(envelope=bbox) %>%
    st_collection_extract()
  #v <- v[unlist(st_intersects(samps, v))]
  v <- st_as_sf(v)
  
  v$Values <- samps$lyr.1
  
  p <- ggplot(v, aes(fill=Values)) +
    geom_sf() +
    scale_fill_binned_divergingx(palette="Geyser") +
    theme_classic() + 
    labs(title=paste(level, "variation habitat"), subtitle="Voronoi tesselation") +
    lims(x=c(-n, n), y=c(-n, n)) +
    guides(fill=guide_legend(title="Value"))
  
  png(paste0("plots/v_plots/v_",level,index,".png"), 6, 6, res=600, units="in")
  print(p)
  dev.off()
  
  #~~~~~~~~~~~~~~~~~~~~
  # Validate method
  #~~~~~~~~~~~~~~~~~~~~
  ### The following code is similar to the grid algorithm validation
  rast.vals <- matrix(values(x), n*2, n*2)
  total.mean <- mean(rast.vals)
  
  ### Find point values from the centroid of each cell
  cents.pts <- crds(x)
  cents <- terra::extract(x, cents.pts)
  cents <- cbind(cents, cents.pts)
  cents <- st_as_sf(cents, coords=c('x', 'y'))
  
  matchy <- st_contains(v, cents)
  touchy <- st_touches(v, cents)
  matchymatchy <- matrix(NA, n*2, n*2)
  
  ssw <- NA
  ssb <- NA
  for(i in 1:n.clust){
    want.these <- c(matchy[[i]], touchy[[i]])
    
    m.want <- mean(rast.vals[want.these])
    
    ssw[i] <- sum((rast.vals[want.these] - m.want)^2)
    ssb[i] <- (m.want - total.mean)^2 * length(want.these)
    
    matchymatchy[want.these] <- i
  }
  ss.total <- sum(ssb + ssw)
  ssb.sst <- sum(ssb) / ss.total
  
  
  new.clusts <- order.clusts(c(matchymatchy))
  if(index == "Dunn"){
    score <- dunn.score(x, new.clusts$new)
  }else{
    score <- unlist(intCriteria(x1, new.clusts$new, crit=index)) 
  }
  
  return(list(new.hab=v, p=p, ss.total=ss.total, ssb.sst=ssb.sst, ssw=ssw, ssb=ssb,
              clusters=c(matchymatchy), score=score))
}