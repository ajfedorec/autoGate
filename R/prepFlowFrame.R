prep.flowFrame <- function(flow_frame) {
  # Remove 0 values as these produce -Inf values when log10 is applied
#   trimmed_flow_frame <- flow_frame
#   for (marker in colnames(flow_frame)) {
#     trimmed_flow_frame <- Subset(trimmed_flow_frame, as.logical(exprs(trimmed_flow_frame[, marker]) > 0.0))
#   }

  trimmed_flow_frame <- flow_frame
  for (marker in c("FSC-H", "SSC-H", "SSC-A")) {
    trimmed_flow_frame <- flowCore::Subset(trimmed_flow_frame, as.logical(flowCore::exprs(trimmed_flow_frame[, marker]) > 0.0))
  }

  # log10 transform the values
  log10_flow_frame <- flowCore::transform(trimmed_flow_frame,
                                          flowCore::transformList(from=flowCore::colnames(trimmed_flow_frame), tfun=log10))
  return(log10_flow_frame)
}
