get.bacteria <- function(flow_frame) {
  ## This adjusts the number of clusters to fit based on if there might be a
  ## cluster for background debris.
  ## If the FACS data has been thresholded below 8000 we see background junk
  ## forming a second cluster.
  ## n.b. the 8000 figure is likely to be changeable, maybe for each run?
  ## n.b. in some data we see 2 clusters of bacteria, so 3 clusters. This
  ##      could cause problems.
  if (min(flowCore::exprs(flow_frame$`FSC-H`)) <= log10(8000)){
    K_clusters <- 2
  }   else {
    K_clusters <- 1; #####CHANGEME: 2 HARDCODED VALUES #####
  }

  ## remove outliers based on forward and side scatter
  ## NOT SURE IF 90% or 95% is better
  ##  This fits a cluster (or 2 clusters) to the data and removes 95% outliers
  res95 <- flowClust::flowClust(flow_frame,
                     varNames = c("FSC-H", "SSC-H"),
                     K = K_clusters,
                     level = 0.95) #####CHANGEME: HARDCODED VALUE #####
  bact_flow_frame <- flowClust::split(flow_frame, res95)[[which.max(res95@mu[, 1])]]

  return(bact_flow_frame)
}
