#' Filter bins based on average per base mappability and GC content (optional).
#'
#' @param heal_list Output list from count_heal_data()
#' @param mappability_threshold Threshold average per bin mappability value below which bins are ignored ('0.9' by default).
#' @param gc_quantile Bins with GC content below first and above last quantiles are ignored. Set to 'FALSE' for no filtering ('FALSE' by default).
#' @param count_threshold How many standard deviations away from the count mean to consider a bin as outlier. Set to 'FALSE' for no count filtering (3 by default).
#'
#' @return A list with one filtered bins data table for each progenitor & the genes data tables if present (any full featureCounts outputs dropped).
#' @export
filter_bins <- function(heal_list, mappability_threshold=0.9, gc_quantile=FALSE, count_threshold=3){

  if(count_threshold!=FALSE){
    smp_means <- get_sample_stats(heal_list,method="mean")
    smp_sd <- get_sample_stats(heal_list,method="sd")

    filtered_list <- lapply(heal_list, function(df){
      current_samples <- setdiff(colnames(df$bins),c("chr", "start", "mappability", "gc_content", "end"))
      for(smp in current_samples){
        threshold_value <- smp_means[[smp]]+count_threshold*smp_sd[[smp]]
        df$bins[[smp]][df$bins[[smp]]>threshold_value] <- NA
      }
        return(list(bins=df$bins,genes=df$genes))
    })
    names(filtered_list) <- names(heal_list)
    heal_list <- filtered_list
  }
  cat("===================\n")
  cat("| Filtering bins:  |\n")
  cat("===================\n")
  cat("\n")
  filtered_list <- foreach::foreach(pr_name=names(heal_list))%do%{

    prog <- heal_list[[pr_name]]

    mappa_keep <- prog$bins$mappability>=mappability_threshold

    if(gc_quantile>0){

      gc_quant <- stats::quantile(prog$bins$gc_content, probs = seq(0, 1, 1/gc_quantile))
      gc_min <- gc_quant[2]
      gc_max <- gc_quant[length(gc_quant)-1]
      gc_keep <- prog$bins$gc_content>=gc_min & prog$bins$gc_content<=gc_max

      which_keep <- gc_keep & mappa_keep

      mappa_perc <- round(sum(mappa_keep)/nrow(prog$bins),2)
      gc_perc <- 1-2*(1/gc_quantile)
      total_perc <- round(sum(which_keep)/nrow(prog$bins),2)
      cat(paste0(pr_name," filtering:  ",total_perc*100,"% passed filters. \n"),
          paste0("                ",gc_perc*100,"% passed ",gc_quantile,"th GC quantile filter. \n"),
          paste0("                ",mappa_perc*100,"% passed ",mappability_threshold," mappability filter. \n"))
      cat("\n")
    } else {

      mappa_perc <- round(sum(mappa_keep)/nrow(prog$bins),2)

      cat(paste0(pr_name," filtering:  ",mappa_perc*100,"% passed ",mappability_threshold," mappability filter. \n"))

      which_keep <- mappa_keep
      cat("\n")
    }

    filtered_df <- prog$bins[which_keep,]

    return(list(bins=filtered_df, genes=prog$genes))

  }

  if(count_threshold!=FALSE){
    cat(paste0("Ignoring bins with counts ",count_threshold," standard deviations above the mean: \n"))

    smp_na <- get_sample_stats(filtered_list,method = "is.na")
    smp_len <- get_sample_stats(filtered_list,method = "length")
    na_freq <- lapply(names(smp_na),function(name){
      round((smp_na[[name]]/smp_len[[name]])*100,2)
    })
    names(na_freq) <- names(smp_na)
    for(i in 1:length(na_freq)){
      smp_name <- names(na_freq)[i]
      cat(paste0("                -",smp_name, " count outlier filtering:  ",na_freq[[smp_name]],"% of remaining bins ignored. \n"))
    }
  }

  names(filtered_list) <- names(heal_list)
  return(filtered_list)

}
