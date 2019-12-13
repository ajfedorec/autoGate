#' Trim .fcs files to remove debris and doublets.
#'
#' \code{trim.fcs} uses mixture models to cluster bacteria from background debris and fits a linear model to SSC-H vs SSC-A to remove doublets.
#'
#' @param dir_path a directory path containing the .fcs files  to be parsed or folders to be recursed through.
#' @param pattern a regex pattern to match particular .fcs files. Default is \code{"*.fcs"} matching all .fcs files.
#' @param flu_channels a list of strings of the fluorescence channels to keep in the trimmed data and plotting. Defaults to "BL1-H".
#' @param calibrate a Boolean flag to determine whether to convert fluorescence to MEF values. Requires an .fcs file with named \code{"*beads*.fcs"}. Defaults to \code{FALSE}.
#' @param do_plot a Boolean flag to determine whether to produce plots showing the trimming of each flowFrame. Defaults to \code{FALSE}.
#' @param pre_cleaned have you pre removed background debris
#' @param cal_bead_peaks a list of lists in the form \code{list(list(channel="BL1-H", peaks=c(0, 200, ...)} of MEF fluorescence values for the calibration beads. Default values for BL1-H and YL2-H.
#'
#' @return nothing is returned. A new folder is created with the trimmed .fcs files and plots if the do_plot flag is TRUE.
#' @export
trim.fcs <- function(dir_path, pattern = "*.fcs", flu_channels=c("BL1-H"),
                     do_plot = F, pre_cleaned = F,
                     calibrate = F, cal_bead_peaks = NA){
  ## Create directory for trimmed flowFrames
  if (!dir.exists(paste(dir_path, "trimmed", sep="_"))) {
    dir.create(paste(dir_path, "trimmed", sep="_"), recursive = T)
  }

  if(calibrate){
    ## First step is to get calibration standard curves
    bead_file <- unlist(list.files(path = dir_path, pattern = utils::glob2rx("*beads*.fcs"),
               full.names = T, recursive = T, include.dirs = T))

    bead_frame <- flowCore::read.FCS(bead_file, emptyValue = F)

    calibration_parameters <- get.calibration(dir_path, bead_frame, flu_channels, cal_bead_peaks)
  }

  all_files <- list.files(path = dir_path, pattern = utils::glob2rx(pattern),
                          full.names = T, recursive = T, include.dirs = T)
  print(paste("Trimming ", length(all_files), " .fcs files.", sep = ""))

  for (next_fcs in all_files) {
    flow_frame <- flowCore::read.FCS(next_fcs, emptyValue = F)

    ## Trim 0 values and log10 transform
    prepped_flow_frame <- prep.flowFrame(flow_frame, flu_channels)

    ## Try to remove background debris by clustering
    bacteria_flow_frame <- get.bacteria(prepped_flow_frame, pre_cleaned)
    try(if(bacteria_flow_frame == 0) {next}, silent = T) # if we haven't found bacteria move on to the next flow frame

    ## Try to remove doublets
    singlet_flow_frame <- get.singlets(bacteria_flow_frame)

    ## Convert to MEF
    if(calibrate){
      out_flow_frame <- to.MEF(singlet_flow_frame, flu_channels, calibration_parameters)
    } else {
      out_flow_frame <- singlet_flow_frame
    }

    ## Save trimmed flowFrames to a new folder
    out_name <- paste(dirname(next_fcs),
                      "_trimmed/",
                      basename(next_fcs),
                      sep = "")
    flowCore::write.FCS(out_flow_frame, out_name)

    ##  Plot a grid of graphs showing the stages of trimming
    if (do_plot) {
      ## Theme
      apatheme = ggplot2::theme_bw() +
        ggplot2::theme(strip.text.x = ggplot2::element_text(size = 8),
                       strip.background = ggplot2::element_rect(colour = "white"),
                       axis.text = ggplot2::element_text(size = 8),
                       axis.text.x = ggplot2::element_text(angle = -40, vjust = 0.5),
                       axis.title = ggplot2::element_text(size = 8),
                       # text = ggplot2::element_text(family = 'Arial'),
                       panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(),
                       panel.border = ggplot2::element_blank(),
                       axis.line = ggplot2::element_line(),
                       legend.title = ggplot2::element_blank())

      ##  This function allows us to take a legend from a plot
      get_legend <- function(myggplot){
        tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(myggplot))
        leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
        legend <- tmp$grobs[[leg]]
        return(legend)
      }

      plts <- list()
      ## FSC-H vs SSC-H
      plt_main <- ggplot2::ggplot() +
        ggplot2::geom_point(data = dplyr::sample_n(as.data.frame(flow_frame[, c("FSC-H", "SSC-H")]@exprs), size = min(2000, flowCore::nrow(flow_frame))),
                            ggplot2::aes(x = log10(`FSC-H`), y = log10(`SSC-H`), color = "all_data"),
                            alpha = 0.1) +
        ggplot2::geom_point(data = dplyr::sample_n(as.data.frame(bacteria_flow_frame[, c("FSC-H", "SSC-H")]@exprs), size = min(2000, flowCore::nrow(bacteria_flow_frame))),
                            ggplot2::aes(x = `FSC-H`, y = `SSC-H`, color = "bacteria"),
                            alpha = 0.1) +
        ggplot2::geom_point(data = dplyr::sample_n(as.data.frame(singlet_flow_frame[, c("FSC-H", "SSC-H")]@exprs), size = min(2000, flowCore::nrow(singlet_flow_frame))),
                            ggplot2::aes(x = `FSC-H`, y = `SSC-H`, color = "single_bacteria"),
                            alpha = 0.1) +
        ggplot2::xlab("log10(FSC-H)") +
        ggplot2::ylab("log10(SSC-H)") +
        ggplot2::xlim(1, 6) +
        ggplot2::ylim(1, 6) + apatheme

      ## Grab the legend to use seperately
      legend <- get_legend(plt_main)
      plt_main <- plt_main + ggplot2::theme(legend.position = "none")

      plts[[1]] <- plt_main

      ## SSC-H vs SSC-A
      plt_single <- ggplot2::ggplot() +
        ggplot2::geom_abline(intercept = 0, slope = 1)+
        ggplot2::geom_point(data = dplyr::sample_n(as.data.frame(flow_frame[, c("SSC-H", "SSC-A")]@exprs), size = min(2000, flowCore::nrow(flow_frame))),
                            ggplot2::aes(x = log10(`SSC-H`), y = log10(`SSC-A`), color = "all_data"),
                            alpha = 0.1) +
        ggplot2::geom_point(data = dplyr::sample_n(as.data.frame(bacteria_flow_frame[, c("SSC-H", "SSC-A")]@exprs), size = min(2000, flowCore::nrow(bacteria_flow_frame))),
                            ggplot2::aes(x = `SSC-H`, y = (`SSC-A`), color = "bacteria"),
                            alpha = 0.1) +
        ggplot2::geom_point(data = dplyr::sample_n(as.data.frame(singlet_flow_frame[, c("SSC-H", "SSC-A")]@exprs), size = min(2000, flowCore::nrow(singlet_flow_frame))),
                            ggplot2::aes(x = `SSC-H`, y = (`SSC-A`), color = "single_bacteria"),
                            alpha = 0.1) +
        ggplot2::xlab("log10(SSC-H)") +
        ggplot2::ylab("log10(SSC-A)") +
        ggplot2::xlim(1, 6) +
        ggplot2::ylim(1, 6)  + apatheme +
        ggplot2::theme(legend.position = "none")

      plts[[2]] <- plt_single

      ## NOTE: the local is necessary here as R is not good at variable scope.
      ## Without the local, each plot gets overwritten by the final plot in the loop.
      ## Also note the "<<-" to assign the plot outside of the local scope.
      for (f.count in 1:length(flu_channels))
        local({
          f.count <- f.count
          filt <- flu_channels[f.count]
          plts[[f.count + 2]] <<- ggplot2::ggplot() +
            ggplot2::geom_density(data = as.data.frame(flow_frame[, filt]@exprs),
                                  ggplot2::aes(x = log10(flow_frame[, filt]@exprs),
                                               y = ..count.., fill = "all_data"),
                                  alpha = 0.5) +
            ggplot2::geom_density(data = as.data.frame(bacteria_flow_frame[, filt]@exprs),
                                  ggplot2::aes(x = bacteria_flow_frame[, filt]@exprs,
                                               y = ..count.., fill = "bacteria"),
                                  alpha = 0.5) +
            ggplot2::geom_density(data = as.data.frame(singlet_flow_frame[, filt]@exprs),
                                  ggplot2::aes(x = singlet_flow_frame[, filt]@exprs,
                                               y = ..count.., fill = "single_bacteria"),
                                  alpha = 0.5) +
            ggplot2::xlab(paste("log10(", filt, ")", sep = "")) +
            ggplot2::xlim(0, 6)  + apatheme +
            ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank()) +
            ggplot2::theme(legend.position = "none")

          if(calibrate){
            plts[[f.count + 2]] <<- plts[[f.count + 2]] +
              ggplot2::geom_density(data = as.data.frame(out_flow_frame[, filt]@exprs),
                                    ggplot2::aes(x = out_flow_frame[, filt]@exprs,
                                                 y = ..count..),
                                    alpha = 0.5, fill = "grey")
          }
        })

      plts[[3 + length(flu_channels)]] <- legend

      plt <- gridExtra::arrangeGrob(grobs = plts, ncol = 2)
      title <- grid::textGrob(paste("Trimming of flow data to remove background and doublets:\n",
                                    flowCore::identifier(flow_frame)), gp = grid::gpar(fontsize = 10))
      padding <- grid::unit(5, "mm")
      plt <- gtable::gtable_add_rows(plt,
                                     heights = grid::grobHeight(title) + padding,
                                     pos = 0)
      plt <- gtable::gtable_add_grob(plt, title, 1, 1, 1, ncol(plt))

      print(paste("Plotting trimmed flowFrame ",
                  flowCore::identifier(flow_frame)))

      ggplot2::ggsave(filename = paste(dirname(out_name),
                                       gsub(".fcs",
                                            "_trimmed.pdf",
                                            basename(out_name)),
                                       sep = "/"),
                      plot = plt,
                      width = 105,
                      height = 50 * ceiling((length(plts)/2)) + 20,
                      units = "mm")
    }
  }
}
