to.MEF <- function(singlet_flow_frame, flu_channels, calibration_parameters){
  out_flow_frame <- singlet_flow_frame
  for(fl in flu_channels){
    fl_idx <- which(calibration_parameters$flu == fl)
    if(length(fl_idx)==1){
      linearTrans <- flowCore::linearTransform(transformationId="Linear-transformation",
                                               a=calibration_parameters[fl_idx,]$m,
                                               b=calibration_parameters[fl_idx,]$b)

      out_flow_frame <- flowCore::transform(out_flow_frame, flowCore::transformList(fl ,linearTrans))
    }
  }

  return(out_flow_frame)
}
