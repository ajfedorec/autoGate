get.calibration <- function(dir_path, bead_frame, flu_channels, cal_bead_peaks){
  if(is.na(cal_bead_peaks)){
    cal_bead_peaks <- list(list(channel="BL1-H", peaks=c(0, 822, 2114, 5911, 17013, 41837, 145365, 287558)),
                            list(channel="YL2-H", peaks=c(0, 218, 581, 1963, 6236, 15267, 68766, 181945)))
  }
  peak_channels <- c()
  for(i in 1:length(cal_bead_peaks)){
    peak_channels <- c(peak_channels, cal_bead_peaks[[i]]$channel)
  }

  ## Log10 transform bead data
  log10_bead_frame <- flowCore::transform(bead_frame,
                                          flowCore::transformList(from = flowCore::colnames(bead_frame), tfun = log10))

  ## gate to remove debris and aggregates
  singlet_cluster <- flowClust::flowClust(log10_bead_frame,
                                        varNames = c("FSC-H", "SSC-H"),
                                        K = 1,
                                        level = 0.8);
  singlet_beads <- log10_bead_frame[singlet_cluster,]

  singlet_plot <- ggplot2::ggplot()+
    ggplot2::geom_point(data=as.data.frame(log10_bead_frame@exprs),
                        ggplot2::aes(`FSC-H`, `SSC-H`), size=0.2, alpha=0.01)+
    ggplot2::geom_density_2d(data=as.data.frame(singlet_beads@exprs), ggplot2::aes(`FSC-H`, `SSC-H`), size=0.1)+
    ggplot2::scale_x_continuous(name = "log10 FSC-H")+
    ggplot2::scale_y_continuous(name = "log10 SSC-H")+
    ggplot2::theme_bw(base_size = 8)
  ggplot2::ggsave(plot = singlet_plot, filename = paste(dir_path,
                        "_trimmed/singlet_beads.pdf",
                        sep = ""), width = 60, height = 60, units = "mm")

  ## cluster on fluorescence channels to find beads
  fluorescence_cluster <- flowClust::flowClust(singlet_beads,
                                               varNames = peak_channels,
                                               randomStart = 10,
                                               K = 8);

  # flowClust::plot(x=fluorescence_cluster, data=singlet_beads)

  calibration_parameters <- c()
  for (i in 1:length(peak_channels)){
    hist_plt <- ggplot2::ggplot()+
      ggplot2::geom_histogram(data=as.data.frame(log10_bead_frame@exprs),
                              ggplot2::aes(log10_bead_frame[,peak_channels[i]]@exprs, y=..count..),
                              bins=200, alpha=0.25)+
      ggplot2::geom_histogram(data=as.data.frame(singlet_beads@exprs),
                              ggplot2::aes(singlet_beads[,peak_channels[i]]@exprs, y=..count..),
                              bins=200, alpha=0.75)+
      ggplot2::geom_vline(xintercept = flowClust::getEstimates(fluorescence_cluster, singlet_beads)$locationsC[,i],
                          linetype=2, size=0.2)+
      ggplot2::scale_x_continuous(paste("log10", peak_channels[i], "(a.u.)"))+
      # ggplot2::scale_y_continuous(trans = "log10")+
      ggplot2::theme_bw(base_size = 8)

    ## get measured peaks and peak mefs (remove first bead due to log10(0))
    peaks <- sort(flowClust::getEstimates(fluorescence_cluster, singlet_beads)$locationsC[,i])[-1]
    cal_peaks <- NA
    for(j in 1:length(cal_bead_peaks)){
      if(cal_bead_peaks[[j]]$channel == peak_channels[i]){
        cal_peaks <- log10(cal_bead_peaks[[j]]$peaks)[-1]
        break
      }
    }
    if(is.na(cal_peaks)){ ## if we don't have calibration data, skip this channel
      next
    }

    ## model: fl_mef = exp(m*log(fl_rfi) + b) - fl_mef_auto
    # model <- nls(cal_peaks ~ exp(m * log(peaks) + b) - fl_mef_auto, start=list(m=1,b=0,fl_mef_auto=0), trace = T)
    model <- nls(cal_peaks ~ m * peaks + b, start=list(m=1,b=0))

    calibration_parameters <- rbind(calibration_parameters,
                                    data.frame(flu=peak_channels[i], m=coef(model)["m"], b=coef(model)["b"]))

    # plot(peaks, cal_peaks)
    new.data <- data.frame(peaks = seq(0,6,len = 100))
    # lines(new.data$peaks,predict(model,newdata = new.data))

    cal_plt <- ggplot2::ggplot()+
      ggplot2::geom_point(ggplot2::aes(x=peaks, y=cal_peaks))+
      ggplot2::geom_line(ggplot2::aes(new.data$peaks,predict(model,newdata = new.data)), size=0.2)+
      ggplot2::scale_x_continuous(name = paste("log10", peak_channels[i], "(a.u.)"), limits = c(0,6.5))+
      ggplot2::scale_y_continuous(name = paste("log10", peak_channels[i], "(MEF)"), limits = c(0,6.5))+
      ggplot2::theme_bw(base_size = 8)

    plt <- gridExtra::arrangeGrob(grobs = list(hist_plt, cal_plt), ncol = 2)
    ggplot2::ggsave(plot = plt, filename = paste(dir_path,
                                                 "_trimmed/", peak_channels[i], "_calibration.pdf",
                                                 sep = ""), width = 130, height = 60, units = "mm")
  }

  return(calibration_parameters)
}
