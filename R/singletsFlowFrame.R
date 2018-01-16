get.singlets <- function(flow_frame){
  ## fit a line to SSC-H and SSC-A
  lm_fit <- stats::lm(flow_frame[, "SSC-H"]@exprs ~ flow_frame[, "SSC-A"]@exprs)

  ## trim points that fall "too far" from the line
  ##      This is an attempt to remove doublets and other cell clumps.
  ##      #####CHANGEME: HARDCODED VALUE #####
  singlet_flow_frame <- flowCore::Subset(flow_frame,
                                         lm_fit$residuals ^ 2 < 0.01)

  return(singlet_flow_frame)
}
