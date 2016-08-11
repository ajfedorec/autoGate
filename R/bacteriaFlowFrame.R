get.bacteria <- function(flow_frame) {
  ## calculate clusters for K=1, K=2 and K=3
  ## 1-3 is pretty arbitrary - could this be changed
  flowClust.res <- flowClust::flowClust(flow_frame,
                             varNames=c("FSC-H", "SSC-H"),
                             K=1:3,
                             criterion="ICL",
                             level = 0.95);

  ## get the results for the K with the best ICL
  flowClust.res <- flowClust.res[[flowClust.res@index]]

  print(paste(flowClust.res@K, "clusters found"))

  ## this is a horrible hack
  ## The bacterial clusters seem to be close to SSC-H == 4 so we find the cluster
  ## which is closest.
  ## clst.indx is the index of the cluster which has its centre closest to ssc-h = 4
  ##    the "box(4, flowClust.res@lambda)" is because a Box-Cox transformation
  ##    needs to be applied to get the value of 4 in the clustering space.
  clst.indx <- which(abs(flowClust.res@mu[, 2]-box(4, flowClust.res@lambda))==min(abs(flowClust.res@mu[, 2]-box(4, flowClust.res@lambda))))
  bact_flow_frame <- flowClust::split(flow_frame, flowClust.res)[[clst.indx]]

  return(bact_flow_frame)
}
