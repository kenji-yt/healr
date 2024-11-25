#' Plotting function for heal_list.
#'
#' @param heal_list Output list in heal format (such as output from count_heal_data())
#' @param quick_view_sample The name of a sample to plot (as character)('FALSE' by default).
#' @param output_dir The name of a directory to write all plots to. Will create one if nonexistent.
#' @param n_cores Number of cores to use ('1' by default).
#' @param prog_ploidy Ploidy of the progenitors (Assumed to be equal. '2' by default)
#' @param plot_cn Logical: plot a line indicating infered copy number ('TRUE' by default in CN has been estimated).
#' @param add_bins Logical: plot counts for each bin ('TRUE' by default; normalized in plot_cn=TRUE).
#' @param colour_map A vector of colours for each progenitor ('c(purple,orange)' by default).
#' @param specific_chr A vector of characters indicating which chromosomes to plot (plots all by default).
#' @param return_list Logical: return a list of plots if quick_view_sample.
#'
#' @return Either nothing or a list of plots.
#' @export
plot_bins <- function(heal_list, quick_view_sample=FALSE, output_dir=FALSE, n_cores=1, prog_ploidy=2, plot_cn=TRUE, add_bins=TRUE, colour_map=c("purple","orange"), specific_chr=FALSE, return_list=FALSE){

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

    if(sum(names(sample_averages)==quick_view_sample)==0){
      cat("ERROR: Sample name not recognized for quick view of alignment. Exiting..")
      return()
    }

    if(output_dir==FALSE){
      cat(paste0("Quickly plotting for ",quick_view_sample,". \n"))
    }else{
      cat(paste0("Saving ", quick_view_sample, "to ", output_dir,"."))
    }
    samples <- quick_view_sample

    smp_type_map <- smp_type_map[smp_type_map$sample==quick_view_sample,]

  }else if(output_dir==FALSE){

    cat("ERROR: no output directory and no 'quick_view_sample' set. One must be set.")
    return()

  }else{

    cat(paste0("Saving all samples and chromosomes to ", output_dir,"."))
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



#' Plotting function for heal alignment
#'
#' @param heal_list
#' @param genespace_dir
#' @param quick_view_sample
#' @param output_dir
#' @param n_cores
#' @param only_anchors
#' @param add_bins Logical: plot counts for each bin ('TRUE' by default).
#' @param prog_ploidy Ploidy of the progenitors (Assumed to be equal. '2' by default)
#' @param colour_map
#' @param specific_chr
#' @param return_list
#'
#' @return
#' @export
#'
#' @examples
plot_alignment <- function(heal_list, alignment, quick_view_sample=FALSE, output_dir=FALSE, n_cores=1, add_bins=TRUE, prog_ploidy=2, colour_map=c("purple","orange"), specific_chr=FALSE, return_list=FALSE){

  polyploid_samples <- names(alignment)

  sample_averages <- get_sample_stats(heal_list)

  progenitors <- names(heal_list)

  names(colour_map) <- progenitors

  if(quick_view_sample!=FALSE){

    if(sum(polyploid_samples==quick_view_sample)==0){
      cat("ERROR: Sample name not recognized (or not polyploid) for quick view of alignment. Exiting..")
      return()
    }

    if(output_dir==FALSE){
      cat(paste0("Quickly plotting for ",quick_view_sample,". \n"))
    }else{
      cat(paste0("Saving ", quick_view_sample, "to ", output_dir,"."))
    }

    polyploid_samples <- quick_view_sample


  }else if(output_dir==FALSE){
    cat("ERROR: no output directory and no 'quick_view_sample' set. One must be set. Exiting..")
    return()
  }else{
    cat(paste0("Plotting all samples and chromosomes to ",output_dir,"."))
  }

  foreach::foreach(smp=polyploid_samples)%do%{

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

        ref_chr_col_name <- paste0("chr_", ref)
        ref_start_col_name <- paste0("start_", ref)
        ref_end_col_name <- paste0("end_", ref)
        ref_cn_col_name <- paste0("cn_", ref)

        which_rows_aln_dt <- alignment[[smp]][[ref_chr_col_name]]==chr

        x <- (alignment[[smp]][[ref_start_col_name]][which_rows_aln_dt] + alignment[[smp]][[ref_end_col_name]][which_rows_aln_dt]) / 2
        y_line <- alignment[[smp]][[ref]][which_rows_aln_dt]
        subgnm_group <- rep(ref, length(x))

        for(alt in alt_gnms){

          alt_cn_col_name <- paste0("cn_", alt)
          y_alt <- alignment[[smp]][[alt]][which_rows_aln_dt]

          subgnm_group <- c(subgnm_group, rep(alt, length(y_alt)))
          x <- c(x, x) # same coordinates
          y_line <- c(y_line, y_alt)
        }

        lines_df <- data.frame(start=x, copy=y_line, subgenome=subgnm_group)

        if(add_bins==TRUE){

          which_rows_bins_dt <- heal_list[[ref]]$bins$chr==chr

          x_pts <- heal_list[[ref]]$bins$start[which_rows_bins_dt]
          y_pts <- heal_list[[ref]]$bins[[smp]][which_rows_bins_dt]
          y_pts <- (y_pts/sample_averages[[smp]])*prog_ploidy

          pts_df <- data.frame(start=x_pts, points=y_pts)

          bin_plot <- ggplot2::ggplot() +
            ggplot2::geom_line(data = lines_df, ggplot2::aes(x = start, y = copy, colour=subgenome), linewidth = 2) +
            ggplot2::geom_point(data = pts_df, ggplot2::aes(x = start, y = points, colour=colour_map[[ref]]), size = 1, alpha = 0.1) +
            ggplot2::theme_minimal() +
            ggplot2::scale_color_manual(values = colour_map) +
            ggplot2::ylim(0, 8)+
            ggplot2::labs(title = paste(chr, ref, smp), x = "Position", y = "Copy Number") +
            ggplot2::theme_bw()

          if(output_dir!=FALSE){
            ggplot2::ggsave(filename=out_file,bin_plot,device="png", width = 6, height = 4, units = "in")
          }else{
            return(bin_plot)
          }

        }else{

          bin_plot <- ggplot2::ggplot() +
            ggplot2::geom_line(data = lines_df, ggplot2::aes(x = start, y = copy, colour=subgenome), linewidth = 2) +
            ggplot2::theme_minimal() +
            ggplot2::scale_color_manual(values = colour_map) +
            ggplot2::ylim(0, 8)+
            ggplot2::labs(title = paste(chr, smp, ref), x = "Position", y = "Copy Number") +
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


plot_pre_db_data <- function(densities, quick_view_sample=FALSE, output_dir=FALSE, show_discordant=FALSE, heal_list=FALSE, alignment=FALSE, ylim_max=FALSE, colour_vec=FALSE, prog_ploidy=2){

  is_align_and_count_data <- is.list(alignment) & is.list(heal_list)
  if(show_discordant==TRUE & !is_align_and_count_data){
    cat("ERROR: show_discordant is TRUE but no valid heal_list or aligment provided. Exiting..")
    return()
  }

  polyploid_samples <- names(densities)
  progenitors <- names(heal_list)

  if(quick_view_sample!=FALSE){

    if(sum(polyploid_samples==quick_view_sample)==0){
      cat("ERROR: Sample name not recognized (or not polyploid) for quick view of alignment. Exiting..")
      return()
    }

    cat(paste0("Quickly plotting for ",quick_view_sample, ". \n"))
    polyploid_samples <- quick_view_sample


  }else if(output_dir==FALSE){
    cat("ERROR: no output directory and no 'quick_view_sample' set. One must be set. Exiting..")
    return()
  }else{
    cat(paste0("Plotting all samples and chromosomes to ",output_dir,"."))
  }

  foreach::foreach(smp=polyploid_samples)%do%{

    if(output_dir!=FALSE){
      smp_dir <- paste0(output_dir,"/",smp,"/")
      dir.create(smp_dir,showWarnings = F)
    }

    if(ylim_max==FALSE){
      if(is.list(heal_list)==FALSE){
        cat("ERROR: No ylim_max and no heal_list. Set at least one. Exiting..")
        return()
      }else{
        counts_vec <- unlist(lapply(heal_list, function(dt){
        return(dt$bins[[smp]])
      }))
      ylim_max <- mean(counts_vec, na.rm = T)+3*sd(counts_vec, na.rm = T)
      }
    }

    cn_labels <- names(densities[[smp]])

    if(colour_vec==FALSE){
      colour_vec <- rainbow(n=length(cn_labels), s = 0.7, v = 1)
    }else{
      if(length(cn_labels)!=length(colour_vec)){
        cat("ERROR: Colour vector length not matching number of copy number categories. Exiting..")
        return()
      }
    }

    names(colour_vec) <- cn_labels

    contour(densities[[smp]][[cn_labels[1]]], ylim=c(-5, ylim_max), col = colour_vec[cn_labels[1]])

    for(i in 2:length(cn_labels)){

      contour(densities[[smp]][[cn_labels[i]]], ylim=c(0, ylim_max), col = colour_vec[cn_labels[i]], add=TRUE)

    }

    if(show_discordant==TRUE){

      cn_col_indx <- grep("cn_", colnames(alignment[[smp]]))
      which_discord <- rowSums(alignment[[smp]][, ..cn_col_indx])!=prog_ploidy*length(heal_list)

      per_prog_dt_list <- foreach::foreach(prog=progenitors)%do%{
        cn_col_name <- paste0("cn_", prog)
        bin_col_name <- paste0("bin_index_", prog)
        col_keep <- c(cn_col_name, bin_col_name)

        cn_bindex_dt_list <- apply(alignment[[smp]][which_discord, ..col_keep], 1, function(row){
          bindex_vec <- as.numeric(unlist(strsplit(row[[bin_col_name]], ",")))
          cn_vec <- rep(row[[cn_col_name]], length(bindex_vec))
          return(data.table::data.table(bin_index=bindex_vec, cn=cn_vec))
        })

        bin_cn_dt <- unique(data.table::rbindlist(cn_bindex_dt_list))

        merge_dt <- data.table::data.table(chr=heal_list[[prog]]$CN$chr[bin_cn_dt$bin_index], start=heal_list[[prog]]$CN$start[bin_cn_dt$bin_index])

        merge_dt <- merge(heal_list[[prog]]$bins, merge_dt)

        col_keep <- c("gc_content", smp)

        return(cbind(merge_dt[, ..col_keep], cn=bin_cn_dt$cn))
    }

    points_dt <- data.table::rbindlist(per_prog_dt_list)

    points_dt$cn <- paste0("cn_", points_dt$cn)

    points(points_dt$gc_content, points_dt[[smp]], col=colour_vec[points_dt$cn], pch=16, cex=1)

  }
}

#plot_heal <- Go check out riparian
