get.calibration <- function(bead_file, flu_channels, MEF_peaks, manual_peaks, bead_dens_bw){
  bead_frame <- flowCore::read.FCS(bead_file, emptyValue = F)

  out_name <- paste(dirname(bead_file),
                    "_trimmed/",
                    basename(unlist(strsplit(bead_file, split = "[.]"))[1]),
                    sep = "")

  ## default peak postions if none provided
  if(is.na(MEF_peaks)){
    MEF_peaks <- list(list(channel="BL1-H", peaks=c(0, 822, 2114, 5911, 17013, 41837, 145365, 287558)),
                      list(channel="YL2-H", peaks=c(0, 218, 581, 1963, 6236, 15267, 68766, 181945)))
  }

  ## make a vector of channels for which to produce calibration parameters
  peak_channels <- c()
  for(i in 1:length(MEF_peaks)){
    if(MEF_peaks[[i]]$channel %in% flu_channels){
      peak_channels <- c(peak_channels, MEF_peaks[[i]]$channel)
    }
  }

  ## Remove events on the boundary that would cause -Inf when log transformed
  bf <- flowCore::boundaryFilter(x=flu_channels, side="lower")
  bounded_bead_filter <- flowCore::filter(bead_frame, bf)
  bounded_bead_frame <- flowCore::Subset(bead_frame, bounded_bead_filter)

  ## Log10 transform bead data
  log10_bead_frame <- flowCore::transform(bounded_bead_frame,
                                          flowCore::transformList(from = flowCore::colnames(bounded_bead_frame), tfun = log10))

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
  ggplot2::ggsave(plot = singlet_plot, filename = paste(out_name,
                                                        "_singlet_beads.pdf",
                                                        sep = ""), width = 60, height = 60, units = "mm")

  calibration_parameters <- c()
  for (i in 1:length(peak_channels)){
    hist_plt <- ggplot2::ggplot()+
      ggplot2::geom_density(data=as.data.frame(log10_bead_frame@exprs),
                            ggplot2::aes(log10_bead_frame[,peak_channels[i]]@exprs),
                            fill="black", alpha=0.25, bw = bead_dens_bw)+
      ggplot2::geom_density(data=as.data.frame(singlet_beads@exprs),
                            ggplot2::aes(singlet_beads[,peak_channels[i]]@exprs),
                            fill="green", alpha=0.75, bw = bead_dens_bw)+
      ggplot2::scale_x_continuous(paste("log10", peak_channels[i], "(a.u.)"))+
      ggplot2::theme_bw(base_size = 8)

    ##### FIND PEAKS #####
    if(is.na(manual_peaks)){
      ## find peaks based on density estimate
      dens_d <- stats::density(singlet_beads@exprs[,peak_channels[i]], bw = bead_dens_bw)
      peak_table <- data.frame(dens_d[c("x", "y")])[c(F, diff(diff(dens_d$y)>=0)<0),] ## get peaks and heights
      peaks <- dplyr::top_n(peak_table, n = length(MEF_peaks[[i]]$peaks), wt = y)$x ## select only

      hist_plt <- hist_plt +
        ggplot2::geom_vline(xintercept = peaks,
                            linetype=2, size=0.2)

      peaks <- peaks[-1]
    } else { ## OR peaks are being manually identified
      peaks <- NA
      for(j in 1:length(manual_peaks)){
        if(manual_peaks[[j]]$channel == peak_channels[i]){
          peaks <- (manual_peaks[[j]]$peaks)[-1]
          break
        }
      }
      if(is.na(peaks)){ ## if we don't have calibration data, skip this channel
        next
      }

      hist_plt <- hist_plt +
        ggplot2::geom_vline(xintercept = peaks, linetype=2, size=0.2)
    }

    cal_peaks <- NA
    for(j in 1:length(MEF_peaks)){
      if(MEF_peaks[[j]]$channel == peak_channels[i]){
        cal_peaks <- log10(MEF_peaks[[j]]$peaks)[-1]
        break
      }
    }
    if(is.na(cal_peaks)){ ## if we don't have calibration data, skip this channel
      next
    }

    ## model: fl_mef = exp(m*log(fl_rfi) + b) - fl_mef_auto
    # model <- nls(cal_peaks ~ exp(m * log(peaks) + b) - fl_mef_auto, start=list(m=1,b=0,fl_mef_auto=0), trace = T)
    model <- stats::nls(cal_peaks ~ m * peaks + b, start=list(m=1,b=0))

    calibration_parameters <- rbind(calibration_parameters,
                                    data.frame(flu=peak_channels[i], m=stats::coef(model)["m"], b=stats::coef(model)["b"]))

    new.data <- data.frame(peaks = seq(0,6,len = 100))

    cal_plt <- ggplot2::ggplot()+
      ggplot2::geom_point(ggplot2::aes(x=peaks, y=cal_peaks))+
      ggplot2::geom_line(ggplot2::aes(new.data$peaks, stats::predict(model, newdata = new.data)), size=0.2)+
      ggplot2::scale_x_continuous(name = paste("log10", peak_channels[i], "(a.u.)"), limits = c(0,6.5))+
      ggplot2::scale_y_continuous(name = paste("log10", peak_channels[i], "(MEF)"), limits = c(0,6.5))+
      ggplot2::theme_bw(base_size = 8)

    plt <- gridExtra::arrangeGrob(grobs = list(hist_plt, cal_plt), ncol = 2)
    ggplot2::ggsave(plot = plt, filename = paste(out_name,
                                                 "_", peak_channels[i], "_calibration.pdf",
                                                 sep = ""), width = 130, height = 60, units = "mm")
  }

  return(calibration_parameters)
}
