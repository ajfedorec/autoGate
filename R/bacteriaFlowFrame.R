get.bacteria <- function(flow_frame, pre_cleaned) {
  if(pre_cleaned){
    return(flow_frame)
  }

  ## calculate clusters for K=1 and K=2
  all_clusters <- flowClust::flowClust(flow_frame,
                             varNames = c("FSC-H", "SSC-H"),
                             K = 1:2,
                             criterion = "ICL",
                             level = 0.90);

  ## get the results for the K with the best ICL
  best.clusters <- all_clusters[[all_clusters@index]]
  print(paste(best.clusters@K, "clusters found"))

  if (best.clusters@K == 2) { ## if there are two cluster
    bact.indx <- which(flowClust::getEstimates(all_clusters[[2]], flow_frame)$locationsC[,2]
                         == max(flowClust::getEstimates(all_clusters[[2]], flow_frame)$locationsC[,2]))

    clst.fsc <- flowClust::getEstimates(all_clusters[[2]], flow_frame)$locationsC[[bact.indx, 1]]
    clst.ssc <- flowClust::getEstimates(all_clusters[[2]], flow_frame)$locationsC[[bact.indx, 2]]

    if ((clst.fsc < 4) && (clst.ssc < 3)) {
      print("only debris found")
      bact_flow_frame <- 0
    } else {
      Ks <- seq(1:2)
      debris.clusts <- setdiff(Ks, bact.indx)
      split_flow_frame <- flowClust::split(flow_frame, best.clusters,
                                           population = list(debris = debris.clusts,
                                                             non.debris = bact.indx))
      bact_flow_frame <- split_flow_frame$non.debris
    }

  }
  else { ## if there is only one cluster
    clst.fsc <- flowClust::getEstimates(best.clusters, flow_frame)$locationsC[[1]]
    clst.ssc <- flowClust::getEstimates(best.clusters, flow_frame)$locationsC[[2]]
    if ((clst.fsc < 4) && (clst.ssc < 3)) {
      print("only debris found")
      bact_flow_frame <- 0
    } else {
      bact_flow_frame <- flow_frame[best.clusters,]
    }
  }

  # ## if there is only one cluster,
  # if (best.clusters@K == 1) {
  #   clst.fsc <- flowClust::getEstimates(best.clusters, flow_frame)$locationsC[[1]]
  #   clst.ssc <- flowClust::getEstimates(best.clusters, flow_frame)$locationsC[[2]]
  #   if ((clst.fsc < 4) && (clst.ssc < 3)) {
  #     print("only debris found")
  #     bact_flow_frame <- 0
  #   } else {
  #     bact_flow_frame <- flow_frame[best.clusters,]
  #   }
  # } else {
  #   debris.indx <- which(flowClust::getEstimates(all_clusters[[2]], flow_frame)$locationsC[,2]
  #                        == min(flowClust::getEstimates(all_clusters[[2]], flow_frame)$locationsC[,2]))
  #
  #   Ks <- seq(1:best.clusters@K)
  #   non.debris.clusts <- setdiff(Ks, debris.indx)
  #   split_flow_frame <- flowClust::split(flow_frame, best.clusters,
  #                                       population = list(debris = debris.indx,
  #                                                         non.debris = non.debris.clusts))
  #   bact_flow_frame <- split_flow_frame$non.debris
  # }

  return(bact_flow_frame)
}
