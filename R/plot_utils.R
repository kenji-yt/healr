name_average_type <- function(cn_list, n_cores=1, method="median", sample_type=FALSE){

  progenitors <- names(cn_list)
  sample_name_list <- lapply(cn_list,function(df){setdiff(colnames(df$bins),c("chr","start","mappability","gc_content","end"))})
  sample_names  <- unlist(sample_name_list)
  polyploids <- names(table(sample_names))[table(sample_names)==2]
  progenitor_samples <- setdiff(sample_names, polyploids)

  samples <- unique(sample_names)

  if(sample_type==TRUE){
    prog_list <- lapply(progenitors,function(pro){
      pro_smpls <- setdiff(sample_name_list[[pro]],polyploids)
      data.frame(sample=pro_smpls ,type=rep(pro,length(pro_smpls)))
    })
    smp_type_df <- do.call(rbind,prog_list)
    smp_type_df<- rbind(smp_type_df,data.frame(sample=polyploids,type=rep("polyploid",length(polyploids))))
    return(smp_type_df)

  }else{
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

    return(sample_averages)
  }
}


plot_bins <- function(cn_list, quick_view_sample=FALSE, output_dir=FALSE, n_cores=1, prog_ploidy=2, plot_cn=TRUE, add_bins=TRUE, colour_map=c("purple","orange")){

  cn_exist <- unlist(lapply(cn_list,function(list){list$CN}))
  if(is.null(cn_exist) && plot_cn==TRUE){
    cat("ERROR: plot_cn set to 'TRUE' but no CN data table found in cn_list. Setting plot_cn to 'FALSE'.")
    plot_cn <- FALSE
  }

  sample_averages <- name_average_type(cn_list)
  progenitors <- names(cn_list)

  names(colour_map) <- progenitors

  smp_type_map <- name_average_type(cn_list,sample_type = TRUE)

  if(quick_view_sample!=FALSE){
    cat(paste0("Quickly plotting for ",quick_view_sample,". Setting output_dir to 'FALSE'."))
    output_dir <- FALSE
    samples <- quick_view_sample
    sample_averages <- sample_averages[quick_view_sample]
    smp_type_map <- smp_type_map[smp_type_map$sample==quick_view_sample,]
  }else if(output_dir==FALSE){
    cat("ERROR: no output directory or no 'quick view sample' set. One must be set.")
    return(NULL)
  }else{
    cat(paste0("Plotting all samples and chromosomes to ",output_dir,"."))
    samples <- names(sample_averages)
  }

  for(smp in samples){

    if(output_dir!=FALSE){
      ref_dir <- paste0(output_dir,"/",smp,"/")
      dir.create(ref_dir,showWarnings = F)
    }

    lapply(cn_list,function(df){colnames(df$bins)
    })

    if(smp_type_map$type[smp_type_map$sample==smp]!="polyploid"){
      current_prog <- smp_type_map$type[smp_type_map$sample==smp]
    }else{
      current_prog <- progenitors
    }
    for(prog in current_prog){

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

        x <- cn_list[[prog]]$bins$start[which_rows]
        y_pts <- cn_list[[prog]]$bins[[smp]][which_rows]

        if(plot_cn==TRUE){
          y_line <- cn_list[[prog]]$CN[[smp]][which_rows]
          if(add_bins==TRUE){
            y_pts <- (y_pts/sample_averages[[smp]])*prog_ploidy
          }

          plot_df <- data.frame(start=x,counts=y_pts,copy=y_line,group=rep(prog,length(y_line)))

          bin_plot <- ggplot2::ggplot() +
            ggplot2::geom_line(data = plot_df, ggplot2::aes(x = x, y = copy, colour=group), linewidth = 2) +
            ggplot2::geom_point(data = plot_df, ggplot2::aes(x = x, y = counts, colour=group), size = 1, alpha = 0.1) +
            ggplot2::theme_minimal() +
            ggplot2::scale_color_manual(values = colour_map) +
            ggplot2::ylim(0, 8)+
            ggplot2::labs(title = paste(chr,smp), x = "Position", y = "Copy Number") +
            ggplot2::theme_bw()

          if(output_dir!=FALSE){
            ggplot2::ggsave(filename=out_file,bin_plot,device="png", width = 6, height = 4, units = "in")
          }else{
            print(bin_plot)
          }

        }else{

          plot_df <- data.frame(start=x,counts=y_pts, group=rep(prog,length(y_pts)))

          bin_plot <- ggplot2::ggplot() +
            ggplot2::geom_point(data = plot_df, ggplot2::aes(x = x, y = counts, colour=group), size = 1, alpha = 1) +
            ggplot2::theme_minimal() +
            ggplot2::scale_color_manual(values = colour_map) +
            ggplot2::labs(title = paste(c,ref), x = "Position", y = "Counts") +
            ggplot2::theme_bw()

          if(output_dir!=FALSE){
            ggplot2::ggsave(filename=out_file,bin_plot,device="png", width = 6, height = 4, units = "in")
          }else{
            print(bin_plot)
          }
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
