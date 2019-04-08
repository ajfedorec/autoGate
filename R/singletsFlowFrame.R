get.singlets <- function(flow_frame){

  ## method using function from openCYto package
  sg <- flowStats::singletGate(flow_frame, area = "SSC-A", height = "SSC-H",
                               prediction_level = 0.95)
  fres <- flowCore::filter(flow_frame, sg)
  singlet_flow_frame <- flowCore::Subset(flow_frame, fres)

  return(singlet_flow_frame)
}
