get.bacteria <- function(flow_frame) {
  ## calculate clusters for K=1, K=2 and K=3
  ## 1-3 is pretty arbitrary - could this be changed
  all_clusters <- flowClust::flowClust(flow_frame,
                             varNames = c("FSC-H", "SSC-H"),
                             K = 1:3,
                             criterion = "ICL",
                             level = 0.95);

  ## get the results for the K with the best ICL
  best.clusters <- all_clusters[[all_clusters@index]]

  print(paste(best.clusters@K, "clusters found"))

  ## this is a horrible hack
  ## For the flow settings I've been using the bacterial clusters seem to be
  ## close to SSC-H == 4 so we find the cluster which is closest.
  ## clst.indx is the index of the cluster which has its centre closest to
  ## ssc-h = 4
  ##    the "box(4, all_clusters@lambda)" is because a Box-Cox transformation
  ##    needs to be applied to get the value of 4 in the clustering space.
  clst.indx <- which(abs(best.clusters@mu[, 2] - flowClust::box(4, best.clusters@lambda))
                     == min(abs(best.clusters@mu[, 2] - flowClust::box(4, best.clusters@lambda))))
  bact_flow_frame <- flowClust::split(flow_frame, best.clusters)[[clst.indx]]

  return(bact_flow_frame)
}
