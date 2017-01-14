#' Trim .fcs files to remove debris and doublets.
#'
#' \code{trim.fcs} uses mixture models to cluster bacteria from background debris and fits a linear model to SSC-H vs SSC-A to remove doublets.
#'
#' @param dir_path a directory path containing the .fcs files  to be parsed or folders to be recursed through.
#' @param pattern a regex pattern to match particular .fcs files. Default is \code{"*.fcs"} matching all .fcs files.
#' @param do_plot a Boolean flag to determine whether to produce plots showing the trimming of each flowFrame. Defaults to \code{FALSE}.
#'
#' @return nothing is returned. A new folder is created with the trimmed .fcs files and plots if the do_plot flag is TRUE.
#' @export
trim.fcs <- function( dir_path, pattern = "*.fcs", do_plot = F ){
  all_files <- list.files(path = dir_path, pattern = utils::glob2rx(pattern),
                          full.names = T, recursive = T, include.dirs = T)
  print(paste("Trimming ", length(all_files), " .fcs files.", sep = ""))

  for (next_fcs in all_files){
    flow_frame <- flowCore::read.FCS(next_fcs, emptyValue = F)

    ## Trim 0 values and log10 transform
    prepped_flow_frame <- prep.flowFrame(flow_frame)

    ## Try to remove background debris by clustering
    bacteria_flow_frame <- get.bacteria(prepped_flow_frame)

    ## Try to remove doublets
    singlet_flow_frame <- get.singlets(bacteria_flow_frame)


    ## Save trimmed flowFrames to a new folder
    out_name <- paste(dirname(next_fcs),
                      "_trimmed/",
                      basename(next_fcs),
                      sep = "")
    if (!dir.exists(dirname(out_name))){
      dir.create(dirname(out_name), recursive = T)
    }
    flowCore::write.FCS(singlet_flow_frame, out_name)

    ##  Plot a grid of graphs showing the stages of trimming
    if (do_plot){
      ##  This function allows us to take a legend from a plot
      get_legend <- function(myggplot){
        tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(myggplot))
        leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
        legend <- tmp$grobs[[leg]]
        return(legend)
      }

      plt_main <- ggplot2::ggplot() +
        ggplot2::geom_point(data = as.data.frame(flow_frame[, c("FSC-H", "SSC-H")]@exprs),
                            ggplot2::aes(x = log10(`FSC-H`), y = log10(`SSC-H`), color = "all_data"),
                            alpha = 0.1) +
        ggplot2::geom_point(data = as.data.frame(bacteria_flow_frame[, c("FSC-H", "SSC-H")]@exprs),
                            ggplot2::aes(x = `FSC-H`, y = `SSC-H`, color = "bacteria"),
                            alpha = 0.1) +
        ggplot2::geom_point(data = as.data.frame(singlet_flow_frame[, c("FSC-H", "SSC-H")]@exprs),
                            ggplot2::aes(x = `FSC-H`, y = `SSC-H`, color = "single_bacteria"),
                            alpha = 0.1) +
        ggplot2::xlim(0, 7) +
        ggplot2::ylim(0, 7)

      plt_fsch <- ggplot2::ggplot() +
        ggplot2::geom_area(data = as.data.frame(flow_frame[, c("FSC-H")]@exprs),
                           ggplot2::aes(x = log10(`FSC-H`), y = ..count.., fill = "all_data"),
                           alpha = 0.5,
                           stat = "bin") +
        ggplot2::geom_area(data = as.data.frame(bacteria_flow_frame[, c("FSC-H")]@exprs),
                           ggplot2::aes(x = `FSC-H`, y = ..count.., fill = "bacteria"),
                           alpha = 0.5,
                           stat = "bin") +
        ggplot2::geom_area(data = as.data.frame(singlet_flow_frame[, c("FSC-H")]@exprs),
                           ggplot2::aes(x = `FSC-H`, y = ..count.., fill = "single_bacteria"),
                           alpha = 0.5,
                           stat = "bin") +
        ggplot2::xlim(0, 7) +
        ggplot2::theme(legend.position = "none")

      plt_ssch <- ggplot2::ggplot() +
        ggplot2::geom_area(data = as.data.frame(flow_frame[, c("SSC-H")]@exprs),
                           ggplot2::aes(x = log10(`SSC-H`), y = ..count.., fill = "all_data"),
                           alpha = 0.5,
                           stat = "bin") +
        ggplot2::geom_area(data = as.data.frame(bacteria_flow_frame[, c("SSC-H")]@exprs),
                           ggplot2::aes(x = `SSC-H`, y = ..count.., fill = "bacteria"),
                           alpha = 0.5,
                           stat = "bin") +
        ggplot2::geom_area(data = as.data.frame(singlet_flow_frame[, c("SSC-H")]@exprs),
                           ggplot2::aes(x = `SSC-H`, y = ..count.., fill = "single_bacteria"),
                           alpha = 0.5,
                           stat = "bin") +
        ggplot2::xlim(0, 7) +
        ggplot2::theme(legend.position = "none")

      if (is.element("FL1-H", flowCore::colnames(flow_frame))){
        plt_bl1h <- ggplot2::ggplot() +
          ggplot2::geom_area(data = as.data.frame(flow_frame[, c("FL1-H")]@exprs),
                             ggplot2::aes(x = log10(`FL1-H`), y = ..count.., fill = "all_data"),
                             alpha = 0.5,
                             stat = "bin") +
          ggplot2::geom_area(data = as.data.frame(bacteria_flow_frame[, c("FL1-H")]@exprs),
                             ggplot2::aes(x = `FL1-H`, y = ..count.., fill = "bacteria"),
                             alpha = 0.5,
                             stat = "bin") +
          ggplot2::geom_area(data = as.data.frame(singlet_flow_frame[, c("FL1-H")]@exprs),
                             ggplot2::aes(x = `FL1-H`, y = ..count.., fill = "single_bacteria"),
                             alpha = 0.5,
                             stat = "bin") +
          ggplot2::xlim(0, 7) +
          ggplot2::theme(legend.position = "none")
      } else {
#         plt_vl1h <- ggplot2::ggplot() +
#           ggplot2::geom_area(data = as.data.frame(flow_frame[, c("VL1-H")]@exprs),
#                              ggplot2::aes(x = log10(`VL1-H`), y = ..count.., fill = "all_data"),
#                              alpha = 0.5,
#                              stat = "bin") +
#           ggplot2::geom_area(data = as.data.frame(bacteria_flow_frame[, c("VL1-H")]@exprs),
#                              ggplot2::aes(x = `VL1-H`, y = ..count.., fill = "bacteria"),
#                              alpha = 0.5,
#                              stat = "bin") +
#           ggplot2::geom_area(data = as.data.frame(singlet_flow_frame[, c("VL1-H")]@exprs),
#                              ggplot2::aes(x = `VL1-H`, y = ..count.., fill = "single_bacteria"),
#                              alpha = 0.5,
#                              stat = "bin") +
#           ggplot2::xlim(0, 7) +
#           ggplot2::theme(legend.position = "none")
#
#         plt_vl2h <- ggplot2::ggplot() +
#           ggplot2::geom_area(data = as.data.frame(flow_frame[, c("VL2-H")]@exprs),
#                              ggplot2::aes(x = log10(`VL2-H`), y = ..count.., fill = "all_data"),
#                              alpha = 0.5,
#                              stat = "bin") +
#           ggplot2::geom_area(data = as.data.frame(bacteria_flow_frame[, c("VL2-H")]@exprs),
#                              ggplot2::aes(x = `VL2-H`, y = ..count.., fill = "bacteria"),
#                              alpha = 0.5,
#                              stat = "bin") +
#           ggplot2::geom_area(data = as.data.frame(singlet_flow_frame[, c("VL2-H")]@exprs),
#                              ggplot2::aes(x = `VL2-H`, y = ..count.., fill = "single_bacteria"),
#                              alpha = 0.5,
#                              stat = "bin") +
#           ggplot2::xlim(0, 7) +
#           ggplot2::theme(legend.position = "none")

        plt_bl1h <- ggplot2::ggplot() +
          ggplot2::geom_area(data = as.data.frame(flow_frame[, c("BL1-H")]@exprs),
                             ggplot2::aes(x = log10(`BL1-H`), y = ..count.., fill = "all_data"),
                             alpha = 0.5,
                             stat = "bin") +
          ggplot2::geom_area(data = as.data.frame(bacteria_flow_frame[, c("BL1-H")]@exprs),
                             ggplot2::aes(x = `BL1-H`, y = ..count.., fill = "bacteria"),
                             alpha = 0.5,
                             stat = "bin") +
          ggplot2::geom_area(data = as.data.frame(singlet_flow_frame[, c("BL1-H")]@exprs),
                             ggplot2::aes(x = `BL1-H`, y = ..count.., fill = "single_bacteria"),
                             alpha = 0.5,
                             stat = "bin") +
          ggplot2::xlim(0, 7) +
          ggplot2::theme(legend.position = "none")

#         plt_bl2h <- ggplot2::ggplot() +
#           ggplot2::geom_area(data = as.data.frame(flow_frame[, c("BL2-H")]@exprs),
#                              ggplot2::aes(x = log10(`BL2-H`), y = ..count.., fill = "all_data"),
#                              alpha = 0.5,
#                              stat = "bin") +
#           ggplot2::geom_area(data = as.data.frame(bacteria_flow_frame[, c("BL2-H")]@exprs),
#                              ggplot2::aes(x = `BL2-H`, y = ..count.., fill = "bacteria"),
#                              alpha = 0.5,
#                              stat = "bin") +
#           ggplot2::geom_area(data = as.data.frame(singlet_flow_frame[, c("BL2-H")]@exprs),
#                              ggplot2::aes(x = `BL2-H`, y = ..count.., fill = "single_bacteria"),
#                              alpha = 0.5,
#                              stat = "bin") +
#           ggplot2::xlim(0, 7) +
#           ggplot2::theme(legend.position = "none")
#
#         plt_yl1h <- ggplot2::ggplot() +
#           ggplot2::geom_area(data = as.data.frame(flow_frame[, c("YL1-H")]@exprs),
#                              ggplot2::aes(x = log10(`YL1-H`), y = ..count.., fill = "all_data"),
#                              alpha = 0.5,
#                              stat = "bin") +
#           ggplot2::geom_area(data = as.data.frame(bacteria_flow_frame[, c("YL1-H")]@exprs),
#                              ggplot2::aes(x = `YL1-H`, y = ..count.., fill = "bacteria"),
#                              alpha = 0.5,
#                              stat = "bin") +
#           ggplot2::geom_area(data = as.data.frame(singlet_flow_frame[, c("YL1-H")]@exprs),
#                              ggplot2::aes(x = `YL1-H`, y = ..count.., fill = "single_bacteria"),
#                              alpha = 0.5,
#                              stat = "bin") +
#           ggplot2::xlim(0, 7) +
#           ggplot2::theme(legend.position = "none")
#
#         plt_yl2h <- ggplot2::ggplot() +
#           ggplot2::geom_area(data = as.data.frame(flow_frame[, c("YL2-H")]@exprs),
#                              ggplot2::aes(x = log10(`YL2-H`), y = ..count.., fill = "all_data"),
#                              alpha = 0.5,
#                              stat = "bin") +
#           ggplot2::geom_area(data = as.data.frame(bacteria_flow_frame[, c("YL2-H")]@exprs),
#                              ggplot2::aes(x = `YL2-H`, y = ..count.., fill = "bacteria"),
#                              alpha = 0.5,
#                              stat = "bin") +
#           ggplot2::geom_area(data = as.data.frame(singlet_flow_frame[, c("YL2-H")]@exprs),
#                              ggplot2::aes(x = `YL2-H`, y = ..count.., fill = "single_bacteria"),
#                              alpha = 0.5,
#                              stat = "bin") +
#           ggplot2::xlim(0, 7) +
#           ggplot2::theme(legend.position = "none")
      }

      plt_single <- ggplot2::ggplot() +
        ggplot2::geom_point(data = as.data.frame(flow_frame[, c("SSC-H", "SSC-A")]@exprs),
                            ggplot2::aes(x = log10(`SSC-H`), y = log10(`SSC-A`), color = "all_data"),
                            alpha = 0.1) +
        ggplot2::geom_point(data = as.data.frame(bacteria_flow_frame[, c("SSC-H", "SSC-A")]@exprs),
                            ggplot2::aes(x = `SSC-H`, y = (`SSC-A`), color = "bacteria"),
                            alpha = 0.1) +
        ggplot2::geom_point(data = as.data.frame(singlet_flow_frame[, c("SSC-H", "SSC-A")]@exprs),
                            ggplot2::aes(x = `SSC-H`, y = (`SSC-A`), color = "single_bacteria"),
                            alpha = 0.1) +
        ggplot2::xlim(0, 7) +
        ggplot2::ylim(0, 7) +
        ggplot2::theme(legend.position = "none")

      legend <- get_legend(plt_main)
      plt_main <- plt_main + ggplot2::theme(legend.position = "none")
#       plt <- gridExtra::arrangeGrob(plt_main, plt_single, legend, plt_ssch,
#                          plt_bl1h, plt_yl1h, ncol = 3, nrow = 2)
#       plt <- gridExtra::arrangeGrob(plt_main, plt_single, legend,
#                                     plt_vl1h, plt_vl2h, plt_bl1h,
#                                     plt_bl2h, plt_yl1h, plt_yl2h,
#                                     ncol = 3, nrow = 3)
        plt <- gridExtra::arrangeGrob(plt_main, plt_single,
                                      plt_bl1h, legend,
                                      ncol = 2, nrow = 2)
      title <- grid::textGrob(paste("Trimming of flow data to remove background
                                    and doublets:\n ",
                                    flowCore::identifier(flow_frame)))
      padding <- grid::unit(5, "mm")
      plt <- gtable::gtable_add_rows(plt,
                                     heights = grid::grobHeight(title) + padding,
                                     pos = 0)
      plt <- gtable::gtable_add_grob(plt, title, 1, 1, 1, ncol(plt))

      ggplot2::ggsave(filename = paste(dirname(out_name),
                              gsub(".fcs",
                                   "_trimmed.png",
                                   basename(out_name)),
                              sep = "/"), plot = plt)
      print(paste("Plotting trimmed flowFrame ",
                  flowCore::identifier(flow_frame)))
    }
  }
}
