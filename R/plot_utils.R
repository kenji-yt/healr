#' Plotting function for heal_list.
#'
#' @param heal_list Output list in heal format (such as output from count_heal_data())
#' @param quick_view_sample The name of a sample to plot (as character)('FALSE' by default).
#' @param output_dir The name of a directory to write all plots to. Will create one if nonexistent.
#' @param n_cores Number of cores to use ('1' by default).
#' @param prog_ploidy Ploidy of the progenitors (Assumed to be equal. '2' by default)
#' @param plot_cn Logical: plot a line indicating infered copy number ('TRUE' by default in CN has been estimated).
#' @param add_bins Logical: plot counts for each bin (normalized in plot_cn=TRUE).
#' @param colour_map A vector of colours for each progenitor ('c(purple,orange)' by default).
#' @param specific_chr A vector of characters indicating which chromosomes to plot (plots all by default).
#' @param return_list Logical: return a list of plots if quick_view_sample.
#'
#' @return Either nothing or a list of plots.
#' @export
plot_bins <- function(heal_list, quick_view_sample=FALSE, output_dir=FALSE, n_cores=1, prog_ploidy=2, plot_cn=TRUE, add_bins=TRUE, colour_map=c("purple","orange"), specific_chr=FALSE, return_list=FALSE){

  cat("Maybe good to have option to save plots for 1 sample only..?")
  cat("To do that just don't set output dir to false If it is provided alongside quick view sample")
  cn_exist <- unlist(lapply(heal_list,function(list){list$CN}))
  if(is.null(cn_exist) && plot_cn==TRUE){
    cat("ERROR: plot_cn set to 'TRUE' but no CN data table found in heal_list. Setting plot_cn to 'FALSE'.")
    plot_cn <- FALSE
  }

  sample_averages <- get_sample_stats(heal_list)
  progenitors <- names(heal_list)

  names(colour_map) <- progenitors

  smp_type_map <- get_sample_stats(heal_list,sample_type = TRUE)

  if(quick_view_sample!=FALSE){
    cat(paste0("Quickly plotting for ",quick_view_sample,". Setting output_dir to 'FALSE'."))
    output_dir <- FALSE
    samples <- quick_view_sample

    smp_type_map <- smp_type_map[smp_type_map$sample==quick_view_sample,]

  }else if(output_dir==FALSE){
    cat("ERROR: no output directory and no 'quick_view_sample' set. One must be set.")
    return()
  }else{
    cat(paste0("Plotting all samples and chromosomes to ",output_dir,"."))
    samples <- names(sample_averages)
  }

  for(smp in samples){

    if(output_dir!=FALSE){
      ref_dir <- paste0(output_dir,"/",smp,"/")
      dir.create(ref_dir,showWarnings = F)
    }

    if(smp_type_map$type[smp_type_map$sample==smp]!="polyploid"){
      current_prog <- smp_type_map$type[smp_type_map$sample==smp]
    }else{
      current_prog <- progenitors
    }

    plot_prog_list <- foreach::foreach(prog=current_prog)%do%{

      if(output_dir!=FALSE){
      dir.create(paste0(ref_dir,"/",prog,"/"),showWarnings = F)
      }

      if(sum(specific_chr!=FALSE)!=FALSE){
        chromo <- intersect(specific_chr,unique(heal_list[[prog]]$bins$chr))
      }else{
        chromo <- unique(heal_list[[prog]]$bins$chr)
      }

      doParallel::registerDoParallel(n_cores)
      plot_list <- foreach::foreach(chr=chromo)%dopar%{

        if(output_dir!=FALSE){
          out_file <- paste0(ref_dir,"/",prog,"/",chr,".png")
        }

        cat("REMOVE THIS LATTER")
        which_bin_rows <- heal_list[[prog]]$bins$chr==chr
        which_cn_rows <- heal_list[[prog]]$CN$chr==chr

        x <- heal_list[[prog]]$bins$start[which_bin_rows]
        y_pts <- heal_list[[prog]]$bins[[smp]][which_bin_rows]

        if(plot_cn==TRUE){
          y_line <- heal_list[[prog]]$CN[[smp]][which_cn_rows]
          if(add_bins==TRUE){
            y_pts <- (y_pts/sample_averages[[smp]])*prog_ploidy
          }

          plot_df <- data.frame(start=x,counts=y_pts,copy=y_line,progenitor=rep(prog,length(y_line)))

          bin_plot <- ggplot2::ggplot() +
            ggplot2::geom_line(data = plot_df, ggplot2::aes(x = x, y = copy, colour=progenitor), linewidth = 2) +
            ggplot2::geom_point(data = plot_df, ggplot2::aes(x = x, y = counts, colour=progenitor), size = 1, alpha = 0.1) +
            ggplot2::theme_minimal() +
            ggplot2::scale_color_manual(values = colour_map) +
            ggplot2::ylim(0, 8)+
            ggplot2::labs(title = paste(chr,smp), x = "Position", y = "Copy Number") +
            ggplot2::theme_bw()

          if(output_dir!=FALSE){
            ggplot2::ggsave(filename=out_file,bin_plot,device="png", width = 6, height = 4, units = "in")
          }else{
            return(bin_plot)
          }

        }else{

          plot_df <- data.frame(start=x,counts=y_pts, progenitor=rep(prog,length(y_pts)))

          bin_plot <- ggplot2::ggplot() +
            ggplot2::geom_point(data = plot_df, ggplot2::aes(x = x, y = counts, colour=progenitor), size = 1, alpha = 1) +
            ggplot2::theme_minimal() +
            ggplot2::scale_color_manual(values = colour_map) +
            ggplot2::labs(title = paste(c,ref), x = "Position", y = "Counts") +
            ggplot2::theme_bw()

          if(output_dir!=FALSE){
            ggplot2::ggsave(filename=out_file,bin_plot,device="png", width = 6, height = 4, units = "in")
          }else{
            return(bin_plot)
          }
        }
      }
      doParallel::stopImplicitCluster()
      if(quick_view_sample!=FALSE){
        if(return_list==TRUE){
          names(plot_list) <- chromo
          return(plot_list)
        }else{
          lapply(plot_list,print)
        }
      }
    }
  }
  if(quick_view_sample!=FALSE){
    if(return_list==TRUE){
      names(plot_prog_list) <- current_prog
      return(plot_prog_list)
    }
  }
}



#' Title
#'
#' @param heal_list
#' @param aln_map
#' @param show_non_anchor
#' @param quick_view_sample
#' @param output_dir
#' @param n_cores
#' @param prog_ploidy
#' @param plot_cn
#' @param add_bins
#' @param colour_map
#' @param specific_chr
#' @param return_list
#'
#' @return
#' @export
#'
#' @examples
plot_alignment <- function(heal_list, aln_map, show_non_anchor=TRUE, quick_view_sample=FALSE, output_dir=FALSE, n_cores=1, prog_ploidy=2, plot_cn=TRUE, add_bins=TRUE, colour_map=c("purple","orange"), specific_chr=FALSE, return_list=FALSE){

  cn_exist <- sum(names(heal_list[[1]])=="CN")!=0
  if(cn_exist!=TRUE){
    cat("ERROR: no CN data. Exiting...")
    return()
  }

  sample_averages <- get_sample_stats(heal_list)
  progenitors <- names(heal_list)

  names(colour_map) <- progenitors

  smp_type_map <- get_sample_stats(heal_list,sample_type = TRUE)

  polyploid_samples <- smp_type_map$sample[smp_type_map$type=="polyploid"]

  if(quick_view_sample!=FALSE){

    if(sum(polyploid_samples==quick_view_sample)==0){
      cat("ERROR: Sample name not recognized (or not polyploid) for quick view of alignment. Exiting..")
      return()
    }

    cat(paste0("Quickly plotting for ",quick_view_sample,". Setting output_dir to 'FALSE'. \n"))
    output_dir <- FALSE
    polyploid_samples <- quick_view_sample


  }else if(output_dir==FALSE){
    cat("ERROR: no output directory and no 'quick_view_sample' set. One must be set.")
    return()
  }else{
    cat(paste0("Plotting all samples and chromosomes to ",output_dir,"."))
  }

  for(smp in polyploid_samples){

    if(output_dir!=FALSE){
      smp_dir <- paste0(output_dir,"/",smp,"/")
      dir.create(smp_dir,showWarnings = F)
    }

    plot_ref_list <- foreach::foreach(ref=progenitors)%do%{

      alt_gnms <- setdiff(progenitors,ref)

      if(output_dir!=FALSE){
        ref_dir <- paste0(smp_dir,"/",ref,"/")
        dir.create(ref_dir,showWarnings = F)
      }


      if(sum(specific_chr!=FALSE)!=FALSE){
        chromo <- intersect(specific_chr,unique(heal_list[[ref]]$bins$chr))
      }else{
        chromo <- unique(heal_list[[ref]]$bins$chr)
      }

      doParallel::registerDoParallel(n_cores)
      plot_list <- foreach::foreach(chr=chromo)%dopar%{

        if(output_dir!=FALSE){
          out_file <- paste0(ref_dir,"/",chr,".png")
        }

        which_rows <- heal_list[[ref]]$bins$chr==chr

        x <- heal_list[[ref]]$bins$start[which_rows]
        y_pts <- heal_list[[ref]]$bins[[smp]][which_rows]
        subgnm_group <- rep(ref, length(x))

        if(plot_cn==TRUE){

          y_line <- heal_list[[ref]]$CN[[smp]][which_rows]

          for(alt in alt_gnms){

            which_alt_row <- aln_map[[ref]][[alt]]$ref_chr==chr

            subgnm_group <- c(subgnm_group, rep(alt, sum(which_alt_row)))

            sub_map <- aln_map[[ref]][[alt]][which_alt_row,]

            x <- c(x, sub_map$ref_bin)

            alt_chrs <- unique(sub_map$alt_chr)
            print("WAWAWA kapav ici ene zfr p fanE")
            #ADFAFAEFQFEWFQFWEFEQ

            for(alt_chr in alt_chrs){

              which_alt_row_map <- sub_map$alt_chr==alt_chr
              x_in_alt <- as.numeric(sub_map[[smp]][which_alt_row_map])

              which_alt_row_CN <- heal_list[[alt]]$CN$chr==alt_chr
              sub_CN_list <- heal_list[[alt]]$CN[which_alt_row_CN,]

              alt_y_line <- c()
              for(alt_x in x_in_alt){
                alt_y_line <- c(alt_y_line, sub_CN_list[[smp]][sub_CN_list$start==alt_x])
              }
              y_line <- c(y_line, alt_y_line)

              if(add_bins==TRUE){

                which_alt_row_count <- heal_list[[alt]]$bins$chr==alt_chr
                sub_count_list <- heal_list[[alt]]$bins[which_alt_row_count,]

                alt_y_pts <- c()
                for(alt_x in x_in_alt){
                  alt_y_pts <- c(alt_y_pts, sub_count_list[[smp]][sub_count_list$start==alt_x])
                }
                y_pts <- c(y_pts, alt_y_pts)
              }
            }
          }

          if(add_bins==TRUE){

            y_pts <- (y_pts/sample_averages[[smp]])*prog_ploidy

            plot_df <- data.frame(start=x,counts=y_pts,copy=y_line,subgenome=subgnm_group)

            bin_plot <- ggplot2::ggplot() +
              ggplot2::geom_line(data = plot_df, ggplot2::aes(x = x, y = copy, colour=subgenome), linewidth = 2) +
              ggplot2::geom_point(data = plot_df, ggplot2::aes(x = x, y = counts, colour=subgenome), size = 1, alpha = 0.1) +
              ggplot2::theme_minimal() +
              ggplot2::scale_color_manual(values = colour_map) +
              ggplot2::ylim(0, 8)+
              ggplot2::labs(title = paste(chr, ref, smp), x = "Position", y = "Copy Number") +
              ggplot2::theme_bw()

            if (show_non_anchor==TRUE) {
              smp_col_name <- paste0("aligned_",smp)
              which_row <- aln_map[[ref]][[alt]]$ref_chr==chr
              which_dtw_count <- aln_map[[ref]][[alt]][[smp_col_name]][which_row]=="dtw_counts"
              dtw_count_zones <- aln_map[[ref]][[alt]]$ref_bin[which_row][which_dtw_count]
              bin_plot <- bin_plot + ggplot2::geom_vline(xintercept = dtw_count_zones, color = "red", alpha=0.05)

              which_dtw_cn <- aln_map[[ref]][[alt]][[smp_col_name]][which_row]=="dtw_cn"
              dtw_cn_zones <- aln_map[[ref]][[alt]]$ref_bin[which_row][which_dtw_cn]
              bin_plot <- bin_plot + ggplot2::geom_vline(xintercept = dtw_cn_zones, color = "purple", alpha=0.05)
            }

            if(output_dir!=FALSE){
              ggplot2::ggsave(filename=out_file,bin_plot,device="png", width = 6, height = 4, units = "in")
            }else{
              return(bin_plot)
            }

        }else{

          plot_df <- data.frame(start=x,copy=y_line,subgenome=subgnm_group)

          bin_plot <- ggplot2::ggplot() +
            ggplot2::geom_line(data = plot_df, ggplot2::aes(x = x, y = copy, colour=subgenome), linewidth = 2) +
            ggplot2::theme_minimal() +
            ggplot2::scale_color_manual(values = colour_map) +
            ggplot2::ylim(0, 8)+
            ggplot2::labs(title = paste(chr,smp), x = "Position", y = "Copy Number") +
            ggplot2::theme_bw()

          if(output_dir!=FALSE){
            ggplot2::ggsave(filename=out_file,bin_plot,device="png", width = 6, height = 4, units = "in")
          }else{
            return(bin_plot)
          }
        }
      }else{

          plot_df <- data.frame(start=x,counts=y_pts, progenitor=rep(prog,length(y_pts)))

          bin_plot <- ggplot2::ggplot() +
            ggplot2::geom_point(data = plot_df, ggplot2::aes(x = x, y = counts, colour=progenitor), size = 1, alpha = 1) +
            ggplot2::theme_minimal() +
            ggplot2::scale_color_manual(values = colour_map) +
            ggplot2::labs(title = paste(c,ref), x = "Position", y = "Counts") +
            ggplot2::theme_bw()

          if(output_dir!=FALSE){
            ggplot2::ggsave(filename=out_file,bin_plot,device="png", width = 6, height = 4, units = "in")
          }else{
            return(bin_plot)
          }
        }
      }
      doParallel::stopImplicitCluster()
      if(quick_view_sample!=FALSE){
        if(return_list==TRUE){
          names(plot_list) <- chromo
          return(plot_list)
        }else{
          lapply(plot_list,print)
        }
      }
    }
  }
  if(quick_view_sample!=FALSE){
    if(return_list==TRUE){
      names(plot_ref_list ) <- progenitors
      return(plot_ref_list )
    }
  }
}




#plot_heal <- Go check out riparian
