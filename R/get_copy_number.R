#' Infer copy number using Circular Binary Segmentation (from DNAcopy)
#'
#' @param counts Output list from count_heal_data() or filter_bins()
#' @param n_cores Number of cores to use ('1' by default)
#' @param prog_ploidy Ploidy of the progenitors (Assumed to be equal. '1' by default)
#' @param method WA CBS p done mean! Method to infered copy number in each segment ('median' or 'mean'. 'median' by default)
#' @param full_output Logical: Do you want to also get the full DNAcopy output ('FALSE' by default)
#'
#' @return  A list with one element per progenitor containing the bins and genes elements of the inpiut plus a data table with infered copy number per bin for each sample and GC and mappability for each bin.
#' @export
#'
#' @importFrom foreach %dopar%
#' @importFrom foreach %do%
get_copy_number <- function(counts, n_cores=1, prog_ploidy=2, method="median", full_output=FALSE){

  cat("wa CBS p done mean...")

  if(intersect(method, c("median", "mean"))==0 || length(method)!=1){
    cat("ERROR: Invalid method input. Choose either 'median' or 'mean'")
    return(NULL)
  }

  progenitors <- names(counts)
  sample_name_per_prog <- lapply(counts,function(df){setdiff(colnames(df$bins),c("chr","start","mappability","gc_content","end"))})
  samples  <- unique(unlist(sample_name_per_prog))

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


  cn_list <- foreach::foreach(pr_name=progenitors)%do%{

    prog <- counts[[pr_name]]$bins

    sample_current_prog <- unlist(sample_name_per_prog[[pr_name]])

    doParallel::registerDoParallel(n_cores)
    sample_CN <- foreach::foreach(smp=sample_current_prog)%dopar%{

      cna.object <- DNAcopy::CNA(chrom = prog$chr, maploc = prog$start, genomdat = prog[[smp]])
      smooth_cna <- DNAcopy::smooth.CNA(cna.object)
      sgmnts <- DNAcopy::segment(smooth_cna)

      estimated_CN <- round((sgmnts$output$seg.mean/sample_averages[[smp]])*prog_ploidy)

      copy_number <- rep(NA,length(prog$chr))

      for(i in 1:nrow(sgmnts$segRows)){

        copy_number[sgmnts$segRows[i,1]:sgmnts$segRows[i,2]] <- rep(estimated_CN[i],length(sgmnts$segRows[i,1]:sgmnts$segRows[i,2]))

      }
      return(copy_number)
    }
    doParallel::stopImplicitCluster()


    cn_df <- data.table::data.table(prog$chr, prog$start, data.frame(sample_CN))
    colnames(cn_df) <- c("chr", "start", sample_current_prog)
    data.table::setkey(cn_df, chr, start)

    if(full_output==TRUE){
      return(list(bins=counts[[pr_name]]$bins, genes=counts[[pr_name]]$genes, CN=cn_df, DNAcopy=sgmnts))
    }else{
      return(list(bins=counts[[pr_name]]$bins, genes=counts[[pr_name]]$genes, CN=cn_df))
    }
  }

  names(cn_list) <- progenitors
  return(cn_list)

}


utils::globalVariables(c("pr_name", "smp"))
