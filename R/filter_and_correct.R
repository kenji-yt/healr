#' Filter bins based on average per base mappability, read count and GC content.
#'
#' @param heal_list List in heal format (such as output from count_heal_data()).
#' @param mappability_threshold Threshold average per bin mappability value below which bins are ignored ('0.9' by default).
#' @param gc_quantile Bins with GC content below first and above last quantiles are ignored. Set to 'FALSE' for no filtering ('FALSE' by default).
#' @param count_threshold How many standard deviations above the count mean to consider a bin as outlier. Set to 'FALSE' for no count filtering (3 by default).
#' @param replace_by_NA Should outliers count values be set to NA ('TRUE' by default). Otherwise outliers are set to the threshold defined by count_threshold.
#'
#' @return A list with one filtered bins data table for each progenitor & the genes data tables if present (any full featureCounts outputs dropped).
#' @export
filter_bins <- function(heal_list, mappability_threshold = 0.9, gc_quantile = FALSE, count_threshold = 3, replace_by_NA = TRUE) {
  
  # Filtering the data
  if (count_threshold != FALSE) {
    smp_means <- get_sample_stats(heal_list, method = "mean")
    smp_sd <- get_sample_stats(heal_list, method = "sd")

    filtered_list <- lapply(heal_list, function(df) {
      current_samples <- setdiff(colnames(df$bins), c("chr", "start", "mappability", "gc_content", "end"))
      for (smp in current_samples) {
        threshold_value <- smp_means[[smp]] + count_threshold * smp_sd[[smp]]
        if (replace_by_NA == TRUE) {
          df$bins[[smp]][df$bins[[smp]] > threshold_value] <- NA
        } else {
          df$bins[[smp]][df$bins[[smp]] > threshold_value] <- threshold_value
        }
      }
      return(list(bins = df$bins))
    })
    names(filtered_list) <- names(heal_list)
    heal_list <- filtered_list
  }
  
  # Creating the output message.
  cat("===================\n")
  cat("| Filtering bins:  |\n")
  cat("===================\n")
  cat("\n")
  filtered_list <- foreach::foreach(pr_name = names(heal_list)) %do% {
    prog <- heal_list[[pr_name]]

    mappa_keep <- prog$bins$mappability >= mappability_threshold

    if (gc_quantile > 0) {
      gc_quant <- stats::quantile(prog$bins$gc_content, probs = seq(0, 1, 1 / gc_quantile))
      gc_min <- gc_quant[2]
      gc_max <- gc_quant[length(gc_quant) - 1]
      gc_keep <- prog$bins$gc_content >= gc_min & prog$bins$gc_content <= gc_max

      which_keep <- gc_keep & mappa_keep

      mappa_perc <- round(sum(mappa_keep) / nrow(prog$bins), 2)
      gc_perc <- 1 - 2 * (1 / gc_quantile)
      total_perc <- round(sum(which_keep) / nrow(prog$bins), 2)
      cat(
        paste0(pr_name, " filtering:  ", total_perc * 100, "% passed filters. \n"),
        paste0("                ", gc_perc * 100, "% passed ", gc_quantile, "th GC quantile filter. \n"),
        paste0("                ", mappa_perc * 100, "% passed ", mappability_threshold, " mappability filter. \n")
      )
      cat("\n")
    } else {
      mappa_perc <- round(sum(mappa_keep) / nrow(prog$bins), 2)

      cat(paste0(pr_name, " filtering:  ", mappa_perc * 100, "% passed ", mappability_threshold, " mappability filter. \n"))

      which_keep <- mappa_keep
      cat("\n")
    }

    filtered_df <- prog$bins[which_keep, ]

    return(list(bins = filtered_df))
  }

  if (count_threshold != FALSE) {
    
    if (replace_by_NA == TRUE) {
      cat(paste0("Replacing bins with counts ", count_threshold, " standard deviations above the mean: \n"))
    }else{
      cat(paste0("Smoothing bins with counts ", count_threshold, " standard deviations above the mean: \n"))
    }
      
    smp_na <- get_sample_stats(filtered_list, method = "is.na")
    smp_len <- get_sample_stats(filtered_list, method = "length")
    na_freq <- lapply(names(smp_na), function(name) {
      round((smp_na[[name]] / smp_len[[name]]) * 100, 2)
    })
    names(na_freq) <- names(smp_na)


    for (i in 1:length(na_freq)) {
      smp_name <- names(na_freq)[i]
      if (replace_by_NA == TRUE) {
        cat(paste0("                -", smp_name, " count outlier filtering:  ", na_freq[[smp_name]], "% of remaining bins set to NA.\n"))
      }else{
        cat(paste0("                -", smp_name, " count outlier filtering:  ", na_freq[[smp_name]], "% of remaining bins set to threshold. \n"))
      }
    }
  }

  names(filtered_list) <- names(heal_list)
  return(filtered_list)
}


#' Correct GC content using loess method
#'
#' @param heal_list List in heal format (such as output from count_heal_data()).
#' @param n_threads Number of threads to use ('1' by default).
#' @param ymax How many standard deviations above the count median should the plot limit be (default = 1).
#' @param point_size Size of points on plots (default = 0.01).
#' @param alpha Transparency of points on plots (default = 0.5).
#' @param linewidth Width of loess regression line (default = 1). 
#'
#' @returns A heal list with the counts in the 'bins' data tables corrected based on GC content effect. 
#' Also returns a list of lots of corrected and uncorrected counts in the 'gc_correction_plots' list item.
#' @export
#'
#' @examples
correct_gc <- function(heal_list, n_threads=1, ymax=1, point_size=0.01, alpha=0.5, linewidth = 1){
  
  
  smp_types <- get_sample_stats(heal_list = heal_list, sample_type = TRUE)
  smp_medians <- get_sample_stats(heal_list = heal_list, method = "median")
  
  polyploid_samples <- smp_types$sample[smp_types$type=="polyploid"]
  
  gc_list <- lapply(heal_list, function(prog_list) {
    return(prog_list$bins$gc_content)
  })
  gc_vec <- unlist(gc_list, use.names = FALSE)
  doParallel::registerDoParallel(n_threads)
  corrected_values_and_plots_list <- foreach::foreach(smp = polyploid_samples)%dopar%{
    
    count_list <- lapply(heal_list, function(prog_list) {
      return(prog_list$bins[[smp]])
    })
    count_vec <- unlist(count_list, use.names = FALSE)
  
    ## Take the 0.01 edge values closest to & including the min and max
    min_round <- round(min(gc_vec), digits = 2)
    max_round <- round(max(gc_vec), digits = 2)
    
    gc_range <- seq(min_round, max_round, 0.01)
    
    min_vec <- c(min(gc_vec)-0.0001, min_round) # the 0.001 is to enable less or equal opperation below
    gc_range[1] <- min_vec_2[which.min(min_vec_2)]
    
    max_vec_2 <- c(max(gc_vec), max_round)
    gc_range[length(gc_range)] <-max_vec_2[which.max(max_vec_2)]
    
    median_list <- foreach::foreach(i = 1:(length(gc_range)-1))%do%{
      return(median(count_vec[gc_vec>gc_range[i] & gc_vec<=gc_range[i+1]]))
    }
    median_vec <- unlist(median_list)
    
    gc_midpoints <- gc_range[1:(length(gc_range)-1)] + 0.005
    model <- loess(median_vec ~ gc_midpoints)
    smooth_counts <- predict(model)
    
    offsets <- smp_medians[[smp]] - smooth_counts 
    
    non_na_ranges <- which(!is.na(median_vec))
    corrected_count_vec <- count_vec
    for(i in 1:length(non_na_ranges)){
      range_index <- non_na_ranges[i]
      bottom <- gc_range[range_index]
      top <- gc_range[range_index+1] 
      corrected_count_vec[gc_vec>bottom & gc_vec<=top] <- count_vec[gc_vec>bottom & gc_vec<=top] + offsets[i]
    }
    
    ggplot_dt <- data.frame(gc_content=gc_vec, read_counts=count_vec)
    gg_line_dt <- data.frame(midpoints=gc_midpoints[!is.na(median_vec)], loess=smooth_counts)
    ggplot_dt_corrected <- data.frame(gc_content=gc_vec, corrected_read_counts=corrected_count_vec)
    
    
    plot_list <- list()
   
    plot_list$raw <- ggplot2::ggplot(ggplot_dt, aes(gc_content, read_counts)) +
      geom_point(size = point_size, alpha = alpha) + 
      geom_line(data = gg_line_dt, aes(x = midpoints, y = loess), color = "blue", linewidth = linewidth) +
      geom_hline(yintercept = smp_medians[[smp]], color = "red", linetype = "dashed", linewidth = 0.5) +
      ylim(0,median(count_vec)+ymax*sd(count_vec)) + 
      labs(
        x = "GC content", 
        y = "Read Counts"
      )
    plot_list$corrected <- ggplot2::ggplot(ggplot_dt_corrected, aes(gc_content, corrected_read_counts)) +
      geom_point(size = point_size, alpha = alpha) + 
      geom_hline(yintercept = smp_medians[[smp]], color = "red", linetype = "dashed", linewidth = 0.5) +
      ylim(0,median(count_vec)+ymax*sd(count_vec)) + 
      labs(
        x = "GC content", 
        y = "Corrected Read Counts"
      )
    
    output_list <- list(corrected_count_vec, plot_list)
    names(output_list) <- c("corrected_counts", "plots")
    return(output_list)
  }
  doParallel::stopImplicitCluster()
    
  names(corrected_values_and_plots_list) <- polyploid_samples
  
  corrected_counts_list <- lapply(corrected_values_and_plots_list, function(smp){
    return(smp$corrected_counts)
  })
  corrected_counts_dt <- data.table::as.data.table(corrected_counts_list)
  colnames(corrected_counts_dt) <- polyploid_samples
  
  dt_nrows <- cumsum(unlist(lapply(heal_list, function(prog_list) {
    return(nrow(prog_list$bins))
  })))
  
  per_prog_dt_list<- foreach::foreach(i = 1:length(dt_nrows))%do%{
    if(i == 1){ 
      start <- 1
    }else{
      start <- dt_nrows[i-1]+1
    }
    end <- dt_nrows[i]
    prog_dt <- corrected_counts_dt[start:end, ]
    return(prog_dt)
  }
  names(per_prog_dt_list) <- names(dt_nrows)
  
  corrected_heal_list <- foreach::foreach(prog =  names(dt_nrows))%do%{
    bins_no_coun_dt <- heal_list[[prog]]$bins[, 1:5]
    corrected_counts_dt <- per_prog_dt_list[[prog]]
    corrected_bins_dt <- cbind(bins_no_coun_dt, corrected_counts_dt)
    
    heal_list[[prog]]$bins <- corrected_bins_dt
    
    return(heal_list[[prog]])
  }
  
  names(corrected_heal_list) <- names(heal_list)
  
  plot_list <- lapply(corrected_values_and_plots_list, function(smp){
    return(smp$plots)
  })
  
  corrected_heal_list$gc_correction_plots <- plot_list
  
  return(corrected_heal_list)
}

