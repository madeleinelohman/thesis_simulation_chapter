
order.clusts <- function(cr){
  new1 <- data.frame(org=as.vector(as.integer(cr)), new=NA)
  clusts <- unique(new1$org)
  for(j in 1:length(clusts)){
    new1$new[which(new1$org == clusts[j])] <- j
  }
  new1 <- new1 %>%
    group_by(new)
  return(new1)
}


