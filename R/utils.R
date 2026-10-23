#' Get general statistics about the each sample.
#'
#' @param heal_list Output list in heal format (such as output from count_heal_data()).
#' @param n_threads Number of threads to use ('1' by default).
#' @param method The type of statistic computed ('median' by default). Options are 'median', 'local_median' (by subgenome), 'mean', 'local_mean', 'sd', 'is.na' and 'length'.
#'
#' @return Information about the samples.
#' @export
#' @importFrom foreach %dopar%
get_sample_stats <- function(heal_list, n_threads = 1, method = "median") {
  
  sample_name_list <- lapply(heal_list, function(df) {
    setdiff(colnames(df$bins), c("chr", "start", "mappability", "gc_content", "end"))
  })
  sample_names <- unlist(sample_name_list)
  
  incomplete <- names(which(table(sample_names) != length(heal_list)))
  if (length(incomplete) > 0) {
    stop(paste("These samples are missing from at least one subgenome:", paste(incomplete, collapse = ", ")))
  }
  
  samples <- unique(sample_names)
  
  doParallel::registerDoParallel(n_threads)
  sample_stats <- foreach::foreach(smp = samples) %dopar% {
    if (method == "median") {
      avg <- stats::median(stats::na.omit(unlist(lapply(heal_list, function(df) {
        df$bins[[smp]]
      }))))
      median_list <- rep(list(avg), length(heal_list))
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
      mean_list <- rep(list(avg), length(heal_list))
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