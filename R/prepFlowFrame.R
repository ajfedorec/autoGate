prep.flowFrame <- function(flow_frame) {
    ## remove samples with 0 values for FSC, SSC and FL1 as these create -Inf
    ## values when log10 transformed
    trimmed_flow_frame <- Subset(flow_frame, as.logical(exprs(flow_frame$`FSC-H`) > 0))
    trimmed_flow_frame <- Subset(trimmed_flow_frame, as.logical(exprs(trimmed_flow_frame$`SSC-H`) > 0))
    trimmed_flow_frame <- Subset(trimmed_flow_frame, as.logical(exprs(trimmed_flow_frame$`FL1-H`) > 0))
    trimmed_flow_frame <- Subset(trimmed_flow_frame, as.logical(exprs(trimmed_flow_frame$`FSC-A`) > 0))
    trimmed_flow_frame <- Subset(trimmed_flow_frame, as.logical(exprs(trimmed_flow_frame$`SSC-A`) > 0))

    ## log10 transform the FSC, SSC, and FL1 readings
    log10_flow_frame <- transform(trimmed_flow_frame,
                                  transformList(c("FSC-H", "SSC-H", "FSC-A",
                                                  "SSC-A", "FL1-H"), log10))
    return(log10_flow_frame)
}