#' Get general statistics about the each sample.
#'
#' @param heal_list Output list in heal format (such as output from count_heal_data()).
#' @param n_threads Number of threads to use ('1' by default).
#' @param method The type of statistic computed ('median' by default). Options are 'median()', 'mean()', 'sd()', 'is.na()' and 'length()'.
#' @param sample_type Logical on choice to return information about wether the sample is polyploid or not.
#'
#' @return Information about the samples.
#' @importFrom foreach %dopar%
get_sample_stats <- function(heal_list, n_threads = 1, method = "median", sample_type = FALSE) {
  
  progenitors <- names(heal_list)
  sample_name_list <- lapply(heal_list, function(df) {
    setdiff(colnames(df$bins), c("chr", "start", "mappability", "gc_content", "end"))
  })
  sample_names <- unlist(sample_name_list)
  polyploids <- names(table(sample_names))[table(sample_names) == 2]
  progenitor_samples <- setdiff(sample_names, polyploids)

  samples <- unique(sample_names)
  
 
  if (sample_type == TRUE) {
    prog_list <- lapply(progenitors, function(pro) {
      pro_smpls <- setdiff(sample_name_list[[pro]], polyploids)
      data.frame(sample = pro_smpls, type = rep(pro, length(pro_smpls)))
    })
    smp_type_df <- do.call(rbind, prog_list)
    smp_type_df <- rbind(smp_type_df, data.frame(sample = polyploids, type = rep("polyploid", length(polyploids))))
    return(smp_type_df)
    
  } else {
    doParallel::registerDoParallel(n_threads)
    sample_stats <- foreach::foreach(smp = samples) %dopar% {
      if (method == "median") {
        avg <- stats::median(stats::na.omit(unlist(lapply(heal_list, function(df) {
          df$bins[[smp]]
        }))))
        median_list <- list(avg, avg)
        names(median_list) <- names(heal_list)
        return(median_list)
      } else if (method == "local_median") {
        avg <- lapply(heal_list, function(df) {
          median(df$bins[[smp]], na.rm = TRUE) 
        })
        return(avg)
      } else if (method == "mean") {
        avg <- mean(stats::na.omit(unlist(lapply(heal_list, function(df) {
          df$bins[[smp]]
        }))))
        mean_list <- list(avg, avg)
        names(mean_list) <- names(heal_list)
        return(mean_list)
      } else if (method == "local_mean") {
        avg <- lapply(heal_list, function(df) {
          mean(df$bins[[smp]], na.rm = TRUE) 
        })
        return(avg)
      } else if (method == "sd") {
        stndi <- stats::sd(stats::na.omit(unlist(lapply(heal_list, function(df) {
          df$bins[[smp]]
        }))))
        return(stndi)
      } else if (method == "is.na") {
        na_count <- sum(is.na(unlist(lapply(heal_list, function(df) {
          df$bins[[smp]]
        }))))
        return(na_count)
      } else if (method == "length") {
        n_bins <- length(unlist(lapply(heal_list, function(df) {
          df$bins[[smp]]
        })))
        return(n_bins)
      } else {
        stop("Unknown stats method. Please input either 'sd', 'mean' or 'median'.")
      }
    }
    doParallel::stopImplicitCluster()

    names(sample_stats) <- samples

    return(sample_stats)
  }
}


get_offset <- function(){
  
  
}
