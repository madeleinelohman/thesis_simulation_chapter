# 
# dunn.score <- function(hab, ordered.clusts){
#   v <- as.polygons(hab, aggregate=F)
#   s <- sf::st_as_sf(v)
# 
#   s$new <- ordered.clusts
#   d <- NA
#   cents <- NA
#   for(p in 1:length(unique(ordered.clusts))){
#     d[p] <- max(units::drop_units(st_distance(s[which(s$new == p),], which="Euclidean")), na.rm=T)
#     cents[p] <- st_centroid(st_union(s[which(s$new == p),]))
#   }
#   cents <- st_distance(st_as_sfc(cents), which="Euclidean")
# 
#   numer <- min(cents[which(cents > 0)])
#   dem <- max(d)
# 
#   return(numer/dem)
# }


dunn.score <- function(hab, ordered.clusts){
  v <- as.polygons(hab, aggregate=F)
  s <- sf::st_as_sf(v)

  s$new <- ordered.clusts

  shoot <- st_centroid(s)
  shoot2 <- as.data.frame(cbind(st_coordinates(shoot), shoot$lyr.1))

  d <- NA
  cents <- NA
  for(p in 1:length(unique(ordered.clusts))){
    want <- which(shoot$new == p)
    d[p] <- max(dist(shoot2[want,]))
    cents[p] <- st_centroid(st_union(s[which(s$new == p),]))
  }

  cents2 <- as.data.frame(st_coordinates(st_as_sfc(cents)))
  cents <- dist(cents2)

  numer <- min(cents[which(cents > 0)])
  dem <- max(d)

  return(numer/dem)
}
