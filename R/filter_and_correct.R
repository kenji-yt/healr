#' Filter bins based on average per base mappability, read count and GC content.
#'
#' @param heal_list List in heal format (such as output from count_heal_data()).
#' @param mappability_threshold Threshold average per bin mappability value below which bins are ignored ('0.9' by default).
#' @param gc_threshold How many standard deviations away from the GC content median to consider a bin as outlier. Set to 'FALSE' for no filtering ('2' by default).
#' @param count_threshold How many standard deviations above the count median to consider a bin as outlier. Set to 'FALSE' for no count filtering ('2' by default).
#' @param replace_by_NA Should outliers count values be set to NA ('TRUE' by default). Otherwise outliers are set to the threshold defined by count_threshold.
#' @param log_file A file to write the filtering statistics to ('FALSE' by default).
#'
#' @return A list with one filtered bins data table for each progenitor & the genes data tables if present (any full featureCounts outputs dropped).
#' @export
filter_bins <- function(heal_list, mappability_threshold = 0.9, gc_threshold = 2,
                        count_threshold = 2, replace_by_NA = TRUE, log_file = FALSE) {
  
  if(log_file!=FALSE){
    sink(log_file, append = TRUE, split = TRUE)
  }
  
  # Filtering the data
  if (count_threshold != FALSE) {
    smp_medians <- get_sample_stats(heal_list, method = "median")
    smp_sd <- get_sample_stats(heal_list, method = "sd")

    filtered_list <- lapply(heal_list, function(df) {
      current_samples <- setdiff(colnames(df$bins), c("chr", "start", "mappability", "gc_content", "end"))
      for (smp in current_samples) {
        threshold_value <- smp_medians[[smp]][[1]] + count_threshold * smp_sd[[smp]]
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
  
  if (gc_threshold != FALSE) {
    gc_vec <- unlist(lapply(heal_list, function(prog_list){
      prog_list$bins$gc_content
    }))
    gc_median <- median(gc_vec)
    gc_sd <- sd(gc_vec)
  }
    
  filtered_list <- foreach::foreach(pr_name = names(heal_list)) %do% {
    prog <- heal_list[[pr_name]]

    mappa_keep <- prog$bins$mappability >= mappability_threshold

    if (gc_threshold != FALSE) {
      gc_min <- gc_median - gc_threshold * gc_sd
      gc_max <- gc_median + gc_threshold * gc_sd
      gc_keep <- prog$bins$gc_content >= gc_min & prog$bins$gc_content <= gc_max

      which_keep <- gc_keep & mappa_keep

      mappa_perc <- round(sum(mappa_keep) / nrow(prog$bins), 2)
      gc_perc <- round(sum(gc_keep) / nrow(prog$bins), 2)
      total_perc <- round(sum(which_keep) / nrow(prog$bins), 2)
      cat(
        paste0(pr_name, " filtering:  ", total_perc * 100, "% passed filters. \n"),
        paste0("                ", gc_perc * 100, "% passed ", gc_threshold, " standard deviation away from median GC filter. \n"),
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

  if(log_file!=FALSE){
    sink()
  }
  
  names(filtered_list) <- names(heal_list)
  return(filtered_list)
  
  
}


#' Correct GC content using loess method
#'
#' @param heal_list List in heal format (such as output from count_heal_data()).
#' @param n_threads Number of threads to use ('1' by default).
#' @param n_windows Number of gc windows in which to compute the medians ('10' by default). 
#' @param loess_span The span variable to pass to the loess() function ('0.75' by default).
#' @param local_normalize Use local medians to normalize ('FALSE' by default). 
#' @param prog_ploidy of the progenitors (Assumed to be equal. '2' by default).
#' @param output_dir Either the name of a directory to write all plots to (will create one if nonexistent) or FALSE (plots are printed only). Defaults to FALSE. 
#' @param ymax The y limit of the plots (default = '4').
#' @param cex Size of points on plots (default = '1').
#' @param pch Point type (default = '1').
#' @param alpha Point transparency (default = '0.3')
#' @param linewidth Width of loess regression line (default = '3'). 
#' @param device Plot device (as argument to ggplot2::ggsave) ('pdf' by default).
#'
#' @returns A heal list with the counts in the 'bins' data tables corrected based on GC content effect. 
#' Also returns a list of lots of corrected and uncorrected counts in the 'gc_correction_plots' list item.
#' @export
#'
#' @examples
correct_gc <- function(heal_list, n_threads=1, n_windows = 10, loess_span = 0.75,
                       prog_ploidy = 2, local_normalize = FALSE, output_dir=FALSE,
                       ymax = 4, cex = 1, pch = 1, alpha = 0.3, linewidth = 3, device = "pdf"){
  
  smp_types <- get_sample_stats(heal_list = heal_list, sample_type = TRUE)
  
  if(local_normalize == TRUE){
    smp_medians <- get_sample_stats(heal_list = heal_list, method = "local_median")
  }else{
    smp_medians <- get_sample_stats(heal_list = heal_list, method = "median")
  }
  
  polyploid_samples <- smp_types$sample[smp_types$type=="polyploid"]
  progenitors <- names(heal_list)
  
  # Make a vector to identify the progenitor
  prog_vec <- unlist(foreach::foreach(i = 1:length(heal_list))%do%{
    return(rep(progenitors[i], nrow(heal_list[[progenitors[i]]]$bins)))
  })
  
  gc_list <- lapply(heal_list, function(prog_list) {
    return(prog_list$bins$gc_content)
  })
  gc_vec <- unlist(gc_list, use.names = FALSE)

  # Normalize and get the relationship over all sample
  doParallel::registerDoParallel(n_threads)
  gc_and_norm_counts_list <- foreach::foreach(smp = polyploid_samples)%dopar%{
    
    norm_counts_vec <- unlist(foreach::foreach(prog = progenitors)%do%{
      normalized_counts <- heal_list[[prog]]$bins[[smp]] / smp_medians[[smp]][[prog]] * prog_ploidy
      return(normalized_counts)
    })
    
    out_table <- data.table::data.table(gc = gc_vec, count = norm_counts_vec, sample = rep(smp, length(gc_vec)), progenitor = prog_vec)
  }
  doParallel::stopImplicitCluster()
  
  gc_and_norm_counts_dt <- data.table::rbindlist(gc_and_norm_counts_list)
  
  # Get regular intervals of GC spaning the range, and get median read count in each interval
  start <- min(gc_and_norm_counts_dt$gc) - 0.000001 # To make more than below
  end <- max(gc_and_norm_counts_dt$gc)
  
  total_range <- end - start
  step <- total_range / (n_windows + 1)
  
  gc_range <- seq(from = start, to = end, by = step)
  
  median_list <- foreach::foreach(i = 1:(length(gc_range)-1))%do%{
    return(median(gc_and_norm_counts_dt$count[gc_vec>gc_range[i] & gc_vec<=gc_range[i+1]], na.rm = TRUE))
  }
  median_vec <- unlist(median_list)
  
  mid_points <- gc_range[1:(length(gc_range)-1)] + (step/2)
  
  loess_model <- loess(median_vec ~ mid_points, span = loess_span)
  smooth_counts <- predict(loess_model)
  
  # Correct counts
  offset <- 2 - smooth_counts
  gc_and_norm_counts_dt$corrected <- gc_and_norm_counts_dt$count
  
  for(i in 1:(length(gc_range)-1)){
    gc_and_norm_counts_dt$corrected[gc_vec>gc_range[i] & gc_vec<=gc_range[i+1]] <- gc_and_norm_counts_dt$corrected[gc_vec>gc_range[i] & gc_vec<=gc_range[i+1]] + offset[i]
  }
  
  # Show plots
  if(output_dir != FALSE){
    if(device == "svg"){
      svg(paste0(output_dir, "/gc_correction.", device), width = 14, height = 7)
    }else if(device == "pdf"){
      pdf(paste0(output_dir, "/gc_correction.", device), width = 14, height = 7)
    }else if(device == "png"){
      png(paste0(output_dir, "/gc_correction.", device), width = 14, height = 7)
    }else{
      stop("Device not valid. Choose from 'pdf', 'png' or 'svg'.")
    }
  }
  
  par(mfrow = c(1, 2))
  plot(gc_and_norm_counts_dt$gc, gc_and_norm_counts_dt$count,
       main = "LOESS Fit of Normalized Read Count ~ GC content", ylim = c(0, ymax),
       xlab = "GC content", ylab = "Normalized Count", pch = pch, col = rgb(0, 0, 0, alpha = alpha), cex = cex)
  abline(h = 2, lty = 2)
  lines(mid_points, smooth_counts, col = "blue", lwd = linewidth)
  segments(
    x0 = mid_points,
    y0 = smooth_counts - (0.03 / 2),  
    x1 = mid_points,
    y1 = smooth_counts + (0.03 / 2),   
    col = "blue", lwd = linewidth + 3)
  
  plot(gc_and_norm_counts_dt$gc, gc_and_norm_counts_dt$corrected,
       main = "GC Corrected Read Count ~ GC content", ylim = c(0, ymax),
       xlab = "GC content", ylab = "Corrected Normalized Count", pch = pch, col = rgb(0, 0, 0, alpha = alpha), cex = cex)
  abline(h = 2, lty = 2)
  
  if(output_dir != FALSE){
    dev.off()
  }
  
  # Use the relationship to normalize counts in all samples
  corrected_bins_dt <- foreach::foreach(prog = progenitors)%do%{
    
    doParallel::registerDoParallel(n_threads)
    corrected_values <- foreach::foreach(smp = polyploid_samples)%dopar%{
      
      norm_corrected <- gc_and_norm_counts_dt$count[gc_and_norm_counts_dt$progenitor == prog & gc_and_norm_counts_dt$sample == smp]
      corrected_vec <- norm_corrected * smp_medians[[smp]][[prog]] / prog_ploidy
      
      return(corrected_vec)
    }
    doParallel::stopImplicitCluster()
    names(corrected_values) <- polyploid_samples
    
    smp_dt <- data.table::as.data.table(corrected_values)
    
    bin_info_dt <- heal_list[[prog]]$bins[, 1:5]
    
    prog_corrected_dt <- cbind(bin_info_dt, smp_dt)
    
    return(prog_corrected_dt)
  }
  names(corrected_bins_dt) <- progenitors
  
  for(prog in progenitors){
    heal_list[[prog]]$bins <- corrected_bins_dt[[prog]]
  }
  
  return(heal_list)
}

#' Replace copy number spans shorter than user defined length by border value.
#'
#' @param heal_list List in heal format (such as output from count_heal_data()).
#' @param max_length The maximal length at which a copy number span is to be replaced.
#' @param n_threads Number of threads to use. The function can process each chromosome in parallele, one subgenome at a time, 
#' such that the maximal sensible value is that of the chromosome number in the subgenome with most chromosomes.
#'
#' @returns A heal list with short CN spans replaced by neighbouring values if these are identical to each other.
#' @export
#'
#' @examples
remove_short_spans <- function(heal_list, max_length = 1, n_threads = 1){
  
  smp_types <- get_sample_stats(heal_list = heal_list, sample_type = TRUE)
  samples <- smp_types$sample
  
  progenitors <- names(heal_list)
  
  out_list <- foreach::foreach(prog = progenitors)%do%{
    
    chromosomes <- unique(heal_list[[prog]]$CN$chr)
    
    corrected_vectors <- foreach::foreach(smp = samples)%do%{
      
      doParallel::registerDoParallel(n_threads)
      corrected_per_chr_vectors <- foreach::foreach(chr = chromosomes)%dopar%{

        rle_smp <- rle(heal_list[[prog]]$CN[[smp]][heal_list[[prog]]$CN$chr==chr])
        
        if(length(rle_smp$lengths)>1){
          
          if(rle_smp$lengths[1] <= max_length && rle_smp$lengths[2] > max_length){
            rle_smp$values[1] <- rle_smp$values[2]}
        
          if(length(rle_smp$lengths)>2){
            
            for (i in 2:(length(rle_smp$lengths) - 1)) {
              if (rle_smp$lengths[i] <= max_length &&
                  rle_smp$lengths[i - 1] > max_length &&
                  rle_smp$lengths[i + 1] > max_length &&
                  rle_smp$values[i - 1] == rle_smp$values[i + 1]) {
                
                rle_smp$values[i] <- rle_smp$values[i - 1]
              }
            }
            if(rle_smp$lengths[length(rle_smp)] <= max_length && rle_smp$lengths[length(rle_smp)-1] > max_length){
                rle_smp$values[length(rle_smp)] <- rle_smp$values[length(rle_smp)-1]
            }
          }
        }
        corrected_vec <- inverse.rle(rle_smp)
        return(corrected_vec)
      }
      doParallel::stopImplicitCluster()
      
      corrected_vec <- unlist(corrected_per_chr_vectors)
      return(corrected_vec)
    }
    names(corrected_vectors) <- samples
    
    # replace
    for(smp in samples){
      heal_list[[prog]]$CN[[smp]] <- corrected_vectors[[smp]]
    }
    return(heal_list[[prog]])
  }
  names(out_list) <- progenitors
  
  return(out_list)
}


