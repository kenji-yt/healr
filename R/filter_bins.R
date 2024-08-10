#' Filter bins based on average per base mappability and GC content (optional).
#'
#' @param counts_list Output list from count_heal_data()
#' @param mappability_threshold Threshold average per bin mappability value below which bins are ignored ('0.9' by default).
#' @param GC_quantile Bins with GC content below first and above last quantiles are ignored. Set to 'FALSE' for no filtering ('FALSE' by default).
#'
#' @return A list with one filtered bins data table for each progenitor & the genes data tables if present (any full featureCounts outputs dropped).
#' @export
filter_bins <- function(counts_list, mappability_threshold=0.9, GC_quantile=FALSE){

  filtered_list <- foreach::foreach(pr_name=names(counts_list))%do%{

    prog <- counts_list[[pr_name]]

    mappa_keep <- prog$bins$mappability>=mappability_threshold

    if(GC_quantile>0){

      gc_quant <- stats::quantile(prog$bins$gc_content, probs = seq(0, 1, 1/GC_quantile))
      gc_min <- gc_quant[2]
      gc_max <- gc_quant[length(gc_quant)-1]
      gc_keep <- prog$bins$gc_content>=gc_min & prog$bins$gc_content<=gc_max

      which_keep <- gc_keep & mappa_keep

      mappa_perc <- round(sum(mappa_keep)/nrow(prog$bins),2)
      gc_perc <- 1-2*(1/GC_quantile)
      total_perc <- round(sum(which_keep)/nrow(prog$bins),2)
      cat(paste0(pr_name," filtering:  ",total_perc,"% passed filters. \n"),
          paste0("                ",gc_perc,"% passed ",GC_quantile,"th GC quantile filter. \n"),
          paste0("                ",mappa_perc,"% passed ",mappability_threshold," mappability filter. \n"))

    } else {

      mappa_perc <- round(sum(mappa_keep)/nrow(prog$bins),2)

      cat(paste0(pr_name," filtering:  ",mappa_perc,"% passed ",mappability_threshold," mappability filter. \n"))

      which_keep <- mappa_keep

    }


    filtered_df <- prog$bins[which_keep,]

    return(list(bins=filtered_df, genes=prog$genes))

  }

  names(filtered_list) <- names(counts_list)
  return(filtered_list)

}
