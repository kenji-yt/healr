get_copy_number <- function(counts, n_cores, prog_ploidy=2, method="median"){


    if(intersect(method, c("median", "mean"))==0 || length(method)!=1){
      cat("ERROR: Invalid method input. Choose either 'median' or 'mean'")
      quit()
    }

    progenitors <- names(counts)
    sample_name_per_prog <- lapply(counts,function(df){setdiff(colnames(df$bins),c("chr","start","mappability","gc_content","end"))})
    samples  <- unique(unlist(sample_name_per_prog))

    doParallel::registerDoParallel(n_cores)
    sample_average <- foreach::foreach(s=samples)%dopar%{
      if(method=="median"){
        avg <- stats::median(na.omit(unlist(lapply(counts, function(df){df$bins[[s]]}))))
        return(avg)
      }else{
        avg <- mean(na.omit(unlist(lapply(counts, function(df){df$bins[[s]]}))))
        return(avg)
      }
    }
    doParallel::stopImplicitCluster()

    names(sample_medians) <- samples

    # for each progenitor
    cn_list <- foreach(p=progenitors)%do%{

      prog <- counts[[p]]$bins

      sample_current_prog <- unlist(sample_name_per_prog[[p]])

      # For each polyploid sample
      registerDoParallel(n_cores)
      sample_CN <- foreach(s=sample_current_prog)%dopar%{

        cna.object <- CNA(chrom = prog$chr,maploc = prog$start,genomdat = prog[[s]])
        smooth_cna <- smooth.CNA(cna.object)
        sgmnts <- segment(smooth_cna)

        estimated_CN <- round((sgmnts$output$seg.mean/sample_medians[[s]])*prog_ploidy)

        copy_number <- rep(NA,length(prog$chr))
        # go through each segment
        for(i in 1:nrow(sgmnts$segRows)){

          copy_number[sgmnts$segRows[i,1]:sgmnts$segRows[i,2]] <- rep(estimated_CN[i],length(sgmnts$segRows[i,1]:sgmnts$segRows[i,2]))

        }
        return(copy_number)
      }
      stopImplicitCluster()  # Stop the parallel back end

      cn_df <- data.table(prog$chr,prog$start,data.frame(sample_CN))
      colnames(cn_df) <- c("chr","start",sample_current_prog)
      setkey(cn_df,chr,start)
      return(list(bins=cn_df,genes=counts[[p]]$genes))
    }
    names(cn_list) <- progenitors
    return(cn_list)
  }

