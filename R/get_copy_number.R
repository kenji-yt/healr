#' Infer copy number using Circular Binary Segmentation (from DNAcopy)
#'
#' @param heal_list Output list in heal format (such as output from count_heal_data() or filter_bins()).
#' @param n_threads Number of threads to use ('1' by default).
#' @param prog_ploidy Ploidy of the progenitors (Assumed to be equal. '2' by default).
#' @param method Which method to use to assign a copy number to each segment ('global', 'local' or 'manual'. 'global' by default). Global uses average over both subgenomes, local uses per subgenome average and manual uses values set by users in 'average_list'. 
#' @param average Average measure used to assign a copy number in each segment ('median' or 'mean'. 'median' by default).
#' @param average_list A list with one element per sample with each containing one element per subgenome with values used to assign copy numbers.
#' @param full_output Logical: Do you want to also get the full DNAcopy output ('FALSE' by default).
#'
#' @return  A list with one element per progenitor containing the bins and genes elements of the inpiut plus a data table with infered copy number per bin for each sample and GC and mappability for each bin.
#' @export
#'
#' @importFrom foreach %dopar%
#' @importFrom foreach %do%
get_copy_number <- function(heal_list, n_threads = 1, prog_ploidy = 2, method = "global", average = "median", average_list = FALSE, full_output = FALSE) {
  
  
  if (intersect(method, c("global", "local", "manual")) == 0 || length(average) != 1) {
    stop("Invalid method input. Choose either 'global', 'local' or 'manual'")
  }
  
  if (intersect(average, c("median", "mean")) == 0 || length(average) != 1) {
    stop("Invalid average input. Choose either 'median' or 'mean'")
  }

  progenitors <- names(heal_list)
  sample_name_per_prog <- lapply(heal_list, function(df) {
    setdiff(colnames(df$bins), c("chr", "start", "mappability", "gc_content", "end"))
  })
  samples <- unique(unlist(sample_name_per_prog))

  if(method == "global"){
    global_averages <- get_sample_stats(heal_list, method = average)
    
  }else if(method == "local"){
    local_averages <- get_sample_stats(heal_list, method = paste0("local_", average))
    
  }else if (method == "manual" & is.list(average_list)){
    if (!identical(names(average_list), samples)){
      stop("Invalid average_list input. When method == 'manual' you must provide a average_list with values of median for each sample subgenome.")
    }
  }else if (normalize == "manual" & is.logical(average_list)){
    stop("No 'average_list' provided. When method == 'manual' you must provide an 'average_list' with values of median for each sample subgenome.")
  }
  

  cn_list <- foreach::foreach(pr_name = progenitors) %do% {
    prog <- heal_list[[pr_name]]$bins

    sample_current_prog <- unlist(sample_name_per_prog[[pr_name]])

    doParallel::registerDoParallel(n_threads)
    sample_CN <- foreach::foreach(smp = sample_current_prog) %dopar% {
      
      cna.object <- DNAcopy::CNA(chrom = prog$chr, maploc = prog$start, genomdat = prog[[smp]])
      smooth_cna <- DNAcopy::smooth.CNA(cna.object)
      # Remove any chromosome with only NA (happens with plastids for example):
      only_na_chromo_list <- foreach::foreach(chromo=unique(smooth_cna$chrom)) %do% {
        return(sum(!is.na(smooth_cna$Sample.1[smooth_cna$chrom==chromo]))==0)
      }
      only_na_chromo <- unique(smooth_cna$chrom)[unlist(only_na_chromo_list)]
     
      smooth_not_na_cna <- smooth_cna
      if(length(only_na_chromo>1)){
        for(rm_chr in only_na_chromo){
          smooth_not_na_cna <- smooth_not_na_cna[smooth_not_na_cna$chrom!=rm_chr,]
        }
      }
      sgmnts <- DNAcopy::segment(smooth_not_na_cna)
      
      if(method == "global"){
        estimated_CN <- round((sgmnts$output$seg.mean / global_averages[[smp]]) * prog_ploidy)
        
      }else if(method == "local"){
        estimated_CN <- round((sgmnts$output$seg.mean / local_averages[[smp]][[pr_name]]) * prog_ploidy)

      }else if(method == "manual"){
        estimated_CN <- round((sgmnts$output$seg.mean / average_list[[smp]][[pr_name]]) * prog_ploidy)
         
      }
      
      copy_number <- rep(NA, length(prog$chr))

      for (i in 1:nrow(sgmnts$segRows)) {
        copy_number[sgmnts$segRows[i, 1]:sgmnts$segRows[i, 2]] <- rep(estimated_CN[i], length(sgmnts$segRows[i, 1]:sgmnts$segRows[i, 2]))
      }
      return(copy_number)
    }
    doParallel::stopImplicitCluster()

    cn_dt <- data.table::data.table(prog$chr, prog$start, prog$end, prog$gc_content, data.frame(sample_CN))
    colnames(cn_dt) <- c("chr", "start", "end", "gc_content", sample_current_prog)
    data.table::setkey(cn_dt, chr, start, end, gc_content)

    if (full_output == TRUE) {
      return(list(bins = heal_list[[pr_name]]$bins, CN = cn_dt, DNAcopy = sgmnts))
    } else {
      return(list(bins = heal_list[[pr_name]]$bins, CN = cn_dt))
    }
  }

  names(cn_list) <- progenitors
  return(cn_list)
}


utils::globalVariables(c("pr_name", "smp"))
