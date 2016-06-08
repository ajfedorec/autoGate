source("R/prepFlowFrame.R")
source("R/bacteriaFlowFrame.R")
source("R/singletsFlowFrame.R")

library(flowClust)
library(flowCore)
library(ggplot2)
library(gridExtra)

#' Trim .fcs files to remove debris and doublets.
#'
#' \code{trim.fcs} uses mixture models to cluster bacteria from background debris and fits a linear model to SSC-H vs SSC-A to remove doublets.
#'
#' @param dir_path a directory path containing the .fcs files  to be parsed or folders to be recursed through.
#' @param pattern a regex pattern to match particular .fcs files. Default is \code{"*.fcs"} matching all .fcs files.
#' @param do_plot a Boolean flag to determine whether to produce plots showing the trimming of each flowFrame. Defaults to \code{FALSE}.
#'
#' @return
#' @export
trim.fcs <- function( dir_path, pattern = "*.fcs", do_plot = F ){
    all_files <- list.files(path = dir_path, pattern = glob2rx(pattern),
                            full.names = T, recursive = T, include.dirs = T)

    for (next_fcs in all_files){
        flow_frame <- read.FCS(next_fcs)

        ## Trim 0 values and log10 transform
        prepped_flow_frame <- prep.flowFrame(flow_frame)

        ## Try to remove background debris by clustering
        bacteria_flow_frame <- get.bacteria(prepped_flow_frame)

        ## Try to remove doublets
        singlet_flow_frame <- get.singlets(bacteria_flow_frame)


        ## Save trimmed flowFrames to a new folder
        out_name <- gsub("//", "_trimmed/", next_fcs)
        if (!dir.exists(dirname(out_name))){
            dir.create(dirname(out_name), recursive = T)
        }
        write.FCS(singlet_flow_frame, out_name)

        ##  Plot a grid of graphs showing the stages of trimming
        if (do_plot){
            ##  This function allows us to take a legend from a plot
            get_legend <- function(myggplot){
                tmp <- ggplot_gtable(ggplot_build(myggplot))
                leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
                legend <- tmp$grobs[[leg]]
                return(legend)
            }

            plt_main <- ggplot() +
                geom_point(data=as.data.frame(flow_frame[,c("FSC-H", "SSC-H")]@exprs), aes(x=log10(`FSC-H`), y=log10(`SSC-H`),color="all_data")) +
                geom_point(data=as.data.frame(bacteria_flow_frame[,c("FSC-H", "SSC-H")]@exprs), aes(x=`FSC-H`, y=`SSC-H`,color="bacteria")) +
                geom_point(data=as.data.frame(singlet_flow_frame[,c("FSC-H", "SSC-H")]@exprs), aes(x=`FSC-H`, y=`SSC-H`,color="single_bacteria"))

            plt_fsch <- ggplot() +
                geom_density(data=as.data.frame(flow_frame[,c("FSC-H", "SSC-H")]@exprs), aes(log10(`FSC-H`),fill="all_data")) +
                geom_density(data=as.data.frame(bacteria_flow_frame[,c("FSC-H", "SSC-H")]@exprs), aes(`FSC-H`,fill="bacteria")) +
                geom_density(data=as.data.frame(singlet_flow_frame[,c("FSC-H", "SSC-H")]@exprs), aes(x=`FSC-H`,fill="single_bacteria"))

            plt_ssch <- ggplot() +
                geom_density(data=as.data.frame(flow_frame[,c("FSC-H", "SSC-H")]@exprs), aes(log10(`SSC-H`),fill="all_data")) +
                geom_density(data=as.data.frame(bacteria_flow_frame[,c("FSC-H", "SSC-H")]@exprs), aes(`SSC-H`,fill="bacteria")) +
                geom_density(data=as.data.frame(singlet_flow_frame[,c("FSC-H", "SSC-H")]@exprs), aes(x=`SSC-H`,fill="single_bacteria"))

            plt_bl1h <- ggplot() +
                geom_density(data=as.data.frame(flow_frame[,c("FL1-H", "SSC-H")]@exprs), aes(log10(`FL1-H`),fill="all_data")) +
                geom_density(data=as.data.frame(bacteria_flow_frame[,c("FL1-H", "SSC-H")]@exprs), aes(`FL1-H`,fill="bacteria")) +
                geom_density(data=as.data.frame(singlet_flow_frame[,c("FL1-H", "SSC-H")]@exprs), aes(x=`FL1-H`,fill="single_bacteria"))

            plt_single <- ggplot() +
                geom_point(data=as.data.frame(flow_frame[,c("SSC-H", "SSC-A")]@exprs), aes(x=log10(`SSC-H`), y=log10(`SSC-A`),color="all_data")) +
                geom_point(data=as.data.frame(bacteria_flow_frame[,c("SSC-H", "SSC-A")]@exprs), aes(x=`SSC-H`, y=(`SSC-A`),color="bacteria")) +
                geom_point(data=as.data.frame(singlet_flow_frame[,c("SSC-H", "SSC-A")]@exprs), aes(x=`SSC-H`, y=(`SSC-A`),color="single_bacteria"))

            legend <- get_legend(plt_main)
            plt_main <- plt_main + theme(legend.position = "none")
            plt_single <- plt_single + theme(legend.position = "none")
            plt_fsch <- plt_fsch + theme(legend.position = "none")
            plt_ssch <- plt_ssch + theme(legend.position = "none")
            plt_bl1h <- plt_bl1h + theme(legend.position = "none")
            plt <- arrangeGrob(plt_main, plt_single, legend, plt_fsch, plt_ssch,
                               plt_bl1h, ncol = 3, nrow = 2)
            title <- textGrob(paste("Trimming of flow data to remove background and doublets:\n ", identifier(flow_frame)))
            padding <- unit(5, "mm")
            plt <- gtable::gtable_add_rows(plt, heights = grobHeight(title) + padding, pos = 0)
            plt <- gtable::gtable_add_grob(plt, title, 1, 1, 1, ncol(plt))

            ggsave(filename = paste(dirname(out_name),
                                    gsub(".fcs",
                                         "_trimmed.png",
                                         basename(out_name)),
                                    sep = "/"), plot = plt)
            print(paste("Plotting trimmed flowFrame ", identifier(flow_frame)))
        }
    }
}
