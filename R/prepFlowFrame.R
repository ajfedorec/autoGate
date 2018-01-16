prep.flowFrame <- function(flow_frame, flus) {

  # Remove channels that aren't being kept
  min_flow_frame <- flow_frame[,c("FSC-H", "SSC-H", "SSC-A", flus)]

  # log10 transform the values
  log10_flow_frame <-
    flowCore::transform(min_flow_frame,
                        flowCore::transformList(from = flowCore::colnames(min_flow_frame), tfun = log10))

  # Remove -Inf values produced when log10 is applied
  trimmed_flow_frame <- log10_flow_frame
  for (marker in c("FSC-H", "SSC-H", "SSC-A", flus)) {
    trimmed_flow_frame <-
      flowCore::Subset(trimmed_flow_frame,
                       is.finite(flowCore::exprs(trimmed_flow_frame[, marker]))[,1])
  }

  return(trimmed_flow_frame)
}
