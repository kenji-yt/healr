plot_bins <- function(cn_list, samples, output_dir=FALSE, n_cores, prog_ploidy=2, plot_cn=TRUE, add_bins=TRUE, colour_map=c("purple","orange"), sample_averages){

  cn_exist <- unlist(lapply(cn_list,function(list){list$CN}))
  if(is.null(cn_exist) && plot_cn==TRUE){
    cat("ERROR: plot_cn set to 'TRUE' but no CN data table found in cn_list. Setting plot_cn to 'FALSE'.")
    plot_cn <- FALSE
  }


  progenitors <- names(cn_list)

  names(colour_map) <- progenitors

  for(smp in samples){

    if(output_dir!=FALSE){
      ref_dir <- paste0(output_dir,"/",smp,"/")
      dir.create(ref_dir,showWarnings = F)
    }

    for(prog in progenitors){

      if(output_dir!=FALSE){
      dir.create(paste0(ref_dir,"/",prog,"/"),showWarnings = F)
      }

      chromo <- unique(cn_list[[prog]]$bins$chr)
      doParallel::registerDoParallel(n_cores)
      foreach::foreach(chr=chromo)%dopar%{

        if(output_dir!=FALSE){
          out_file <- paste0(ref_dir,"/",prog,"/",chr,".png")
        }

        which_rows <- cn_list[[prog]]$bins$chr==chr

        x <- cn_list[[prog]]$bins$start
        y_pts <- cn_list[[prog]]$bins[[smp]]

        if(plot_cn==TRUE){
          y_line <- cn_list[[prog]]$CN[[smp]]
          if(add_bins==TRUE){
            y_pts <- y_pts/sample_averages[smp]*prog_ploidy
          }

          plot_df <- data.frame(start=x,counts=y_pts,copy=y_line)

          plot <- ggplot() +
            geom_line(data = plot_df, aes(x = x, y = copy, color = colour_map[prog]), size = 2) +
            geom_point(data = plot_df, aes(x = x, y = counts, color = colour_map[prog]), size = 1, alpha = 0.1) +
            theme_minimal() +
            labs(title = paste(c,ref), x = "Position", y = "Copy Number") +
            theme_bw()

          if(output_dir!=FALSE){
            ggsave(filename=out_file,plot,device="png", width = 6, height = 4, units = "in")
          }
          plot

        }else{

          plot_df <- data.frame(start=x,counts=y_pts)

          plot <- ggplot() +
            geom_point(data = plot_df, aes(x = x, y = counts, color = colour_map[prog]), size = 1, alpha = 1) +
            theme_minimal() +
            labs(title = paste(c,ref), x = "Position", y = "Counts") +
            theme_bw()

          if(output_dir!=FALSE){
            ggsave(filename=out_file,plot,device="png", width = 6, height = 4, units = "in")
          }
          plot
        }
      }
      doParallel::stopImplicitCluster()
    }
  }
}


plot_alignment <- function(cn_list, aln_map, polyploids, output_dir=FALSE, n_cores, prog_ploidy=2, bin_size=NULL, plot_cn=TRUE, add_bins=TRUE, colour_map=c("purple","orange"), sample_averages){

  cn_exist <- unlist(lapply(cn_list,function(list){list$CN}))
  if(is.null(cn_exist) && plot_cn==TRUE){
    cat("ERROR: plot_cn set to 'TRUE' but no CN data table found in cn_list. Setting plot_cn to 'FALSE'.")
    plot_cn <- FALSE
  }else if(is.null(bin_size) && plot_cn==TRUE && add_bins==TRUE){
    cat("ERROR: plot_cn set to 'TRUE' but no bin size is defined. Setting add_bins set to 'FALSE'.")
    add_bins <- FALSE
  }

  progenitors <- names(cn_list)

  names(colour_map) <- progenitors

  for(ref in progenitors){

    alt <- setdiff(progenitors,ref)

    ref_count <- cn_list[[ref]]$bins

    alt_merged_count <- merge(cn_list[[alt]]$bins,aln_map[[alt]],by=c("chr","start"))
  }
} #### UNTIL HERE kind of done
#
#     if(plot_cn==TRUE){
#       ref_cn <- cn_list[[ref]]$CN
#       alt_merged_cn <- merge(aln_map[[alt]],CN_CBS[[alt]],by=c("chr","start"))
#     }
#
#     # Go through each sample
#     for(ply in polyploids){
#
#       if(output_dir!=FALSE){
#         ref_dir <- paste0(output_dir,"/",ply,"/",ref)
#         dir.create(ref_dir,showWarnings = F)
#       }
#
#       chromo <- unique(ref_cn$chr)
#       doParallel::registerDoParallel(n_cores)
#       foreach(chr=chromo)%dopar%{
#
#         if(output_dir!=FALSE){
#           out_file <- paste0(ref_dir,"/",chr,".png")
#         }
#
#         x <- c(ref_count$start[ref_count$chr==chr]+(bin_size/2),alt_merged_count$map_pos[alt_merged_count$map_chr==chr])
#         y_pts <- c(ref_count[[ply]][ref_count$chr==chr],alt_merged_count[[ply]][alt_merged_count$map_chr==chr])
#
#         if(plot_cn==TRUE){
#           y_pts <- y_pts/sample_medians[[ply]]*prog_ploidy # convert to CN
#           y_line <- c(ref_cn[[ply]][ref_cn$chr==chr],alt_merged_cn[[ply]][alt_merged_cn$map_chr==chr])
#         }
#
#         prog_group <- c(rep(ref,sum(ref_count$chr==chr)),rep(alt,sum(alt_merged_count$map_chr==chr)))
#
#         plot_df <- data.frame(start=x,counts=y_pts,copy=y_line,group=prog_group)
#
#         colour_map <- c("purple","orange")
#         names(colour_map) <- c(ref,alt)
#
#         # Plot using ggplot2
#         plot <- ggplot() +
#           # Lines: define data and aes separately for lines
#           geom_line(data = plot_df, aes(x = x, y = copy, group = group, color = group), size = 2) +
#           # Points: define data and aes separately for points
#           geom_point(data = plot_df, aes(x = x, y = counts, color = group), size = 1, alpha = 0.1) +
#           # Manual color scale
#           scale_color_manual(values = colour_map) +
#           # Theme and labels
#           theme_minimal() +
#           labs(title = paste(c,ref), x = "Position", y = "Copy Number") +
#           theme(legend.position = "bottom")+
#           theme_bw()
#
#         # Save the plot
#         ggsave(filename=out_file,plot,device="png", width = 6, height = 4, units = "in")
#       }
#       CLOSE BACKEND
#     }
#   }
# }

plot_heal <- function(cn_list, aln_map=FALSE, output_dir=FALSE, n_cores, prog_ploidy=2, bin_size){

  progenitors <- names(cn_list)
  sample_names  <- unlist(lapply(cn_list,function(df){setdiff(colnames(df$bins),c("chr","start","mappability","gc_content","end"))}))
  polyploids <- names(table(sample_names))[table(sample_names)==2]
  progenitors <- setdiff(sample_names, polyploids)
  samples <- unique(sample_names)

  doParallel::registerDoParallel(n_cores)
  sample_averages <- foreach::foreach(smp=samples)%dopar%{
    if(method=="median"){
      avg <- stats::median(stats::na.omit(unlist(lapply(counts, function(df){df$bins[[smp]]}))))
      return(avg)
    }else{
      avg <- mean(stats::na.omit(unlist(lapply(counts, function(df){df$bins[[smp]]}))))
      return(avg)
    }
  }
  doParallel::stopImplicitCluster()

  names(sample_averages) <- samples

  if(output_dir!=FALSE){
    plots_dir <- paste0(output_dir,"/ratio_plots")
    dir.create(plots_dir,showWarnings = F)
    for(smp in samples){dir.create(paste0(plots_dir,"/",smp),showWarnings = F)}
  }

}
