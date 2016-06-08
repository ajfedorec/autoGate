prep.flowFrame <- function(flow_frame) {

  # Remove 0 values as these produce -Inf values when log10 is applied
  trimmed_flow_frame <- flow_frame
  for (marker in colnames(flow_frame)) {
    trimmed_flow_frame <- Subset(trimmed_flow_frame, as.logical(exprs(trimmed_flow_frame[, marker]) > 0))
  }

  # log10 transform the values
  log10_flow_frame <- transform(trimmed_flow_frame,
                                transformList(colnames(trimmed_flow_frame), log10))
  return(log10_flow_frame)
}
