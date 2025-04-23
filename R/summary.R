#' Summary of the copy number analysis
#'
#' @param heal_list A heal alignment object created with get_heal_alignment().
#' @param n_threads Number of threads to use ('1' by default).
#'
#' @return A list with one entry per sample. Each sample in turn contains a list with information about per progenitor and per chromosome. 
#' The information at each chromosome is run length encoded information of copy number and the distribution of bins in each copy number. 
#' There are also total distributions for each progenitor and for the sample overall.  
#' @export
#'
summarize_cn <- function(heal_list, n_threads=1){
  
  cn_exist <- sum(names(heal_list[[1]]) == "CN") != 0
  if (cn_exist != TRUE) {
    stop("No CN data. Exiting...")
  }
  
  progenitors <- names(heal_list)
  samples <- names(get_sample_stats(heal_list))
  
  summary_per_sample_list <- foreach::foreach(smp=samples)%do%{
    
    which_prog <- unlist(lapply(heal_list, function(prog_list){
      return(sum(colnames(prog_list$CN)==smp)==1)
    }))
    
    valid_progenitors <- progenitors[which_prog]
    
    per_prog_dt_list <- foreach::foreach(prog=valid_progenitors)%do%{
      
      chromo <- unique(heal_list[[prog]]$CN$chr)
      
      doParallel::registerDoParallel(n_threads)
      summary_list <- foreach::foreach(chr=chromo)%dopar%{
        which_row <- heal_list[[prog]]$CN$chr==chr
        counts <- table(heal_list[[prog]]$CN[[smp]][which_row])
        percentage <- counts/sum(counts)*100
        count_table <- rbind(counts, percentage)
        run_length_encoding <- rle(heal_list[[prog]]$CN[[smp]][which_row])
        output_list <- list(count_table=count_table, run_length_encoding=run_length_encoding)
        return(output_list)
      }
      doParallel::stopImplicitCluster()
      names(summary_list) <- chromo
      
      counts <- table(heal_list[[prog]]$CN[[smp]])
      percentage <- counts/sum(counts)*100
      count_table <- rbind(counts, percentage)
      summary_list[[paste0("total_count_table_", prog)]] <- count_table
      return(summary_list)
    }
    names(per_prog_dt_list) <- valid_progenitors
    
    totals <- lapply(per_prog_dt_list, function(prog){
      which_total <- grep("total", names(prog))
      return(prog[[which_total]])
    })
    unique_CNs <- unique(unlist(lapply(totals, colnames)))
    counts_by_cn <- foreach::foreach(cn=unique_CNs)%do%{
      count_list <- lapply(totals, function(dt){
        which_col <- colnames(dt)==cn
        counts <- dt["counts", which_col]
        return(counts)
      })
      total_counts <- sum(unlist(count_list))
    }
    names(counts_by_cn) <- unique_CNs
    counts <- unlist(counts_by_cn)
    percentage <- counts/sum(counts)*100
    count_table <- rbind(counts, percentage)
    
    per_prog_dt_list[["total_count_table"]] <- count_table
    
    return(per_prog_dt_list)
  }
  names(summary_per_sample_list) <- samples
  return(summary_per_sample_list)
}

#' Summary of the alignment outcome for each polyploid sample
#'
#' @param alignment A heal alignment object created with get_heal_alignment().
#' @param n_threads Number of threads to use ('1' by default).
#'
#' @return A list with one entry per polyploid sample. Each sample in turn contains a list with information about per progenitor and per chromosome.
#' The values are both a data table giving run length encoded style information along the chromosomes and another data table giving the proportion of the chromosomes 
#' in each combination of copy number. There are also total distributions for each progenitor and for the sample overall. 
#' @export
#'
summarize_aln <- function(alignment, n_threads=1){
  
  polyploid_samples <- names(alignment)
  
  alignment_summary_list <- foreach::foreach(smp=polyploid_samples)%do%{
    
    chr_col_names <- grep("chr", colnames(alignment[[smp]]), value = TRUE)
    progenitors <- gsub("chr_","", chr_col_names)
    
    names(chr_col_names) <- progenitors
    
    cn_col_names <- gsub("chr", "cn", chr_col_names)
    ratio_name <- paste(names(cn_col_names), collapse="<:>")
    
    summary_per_prog_list <- foreach::foreach(prog=progenitors)%do%{
      
      chromo <- unique(alignment[[smp]][[chr_col_names[[prog]]]])
      doParallel::registerDoParallel(n_threads)
      per_chromo_list <- foreach::foreach(chr=chromo)%dopar%{
        
        which_row <- alignment[[smp]][[chr_col_names[[prog]]]]==chr
        
        start_end_ref_colnames <- paste0(c("start_", "end_"), prog)
        sub_aln_dt <- alignment[[smp]][which_row,]
        sub_aln_dt <- sub_aln_dt[order(sub_aln_dt[[start_end_ref_colnames[1]]]),]
        
        sub_aln_dt[, combination := do.call(paste, c(.SD[, cn_col_names, with = FALSE], sep = "<:>"))]
        
        rle_combination <- rle(sub_aln_dt$combination)
        
        end_indexes <- cumsum(rle_combination$lengths)
        start_regime_vec <- min(sub_aln_dt[1:end_indexes[1], get(start_end_ref_colnames[1])])
        end_regime_vec <- max(sub_aln_dt[1:end_indexes[1], get(start_end_ref_colnames[2])])
        
        for(i in 2:length(end_indexes)){
          
          start_regime <- min(sub_aln_dt[(end_indexes[(i-1)]+1):end_indexes[i], get(start_end_ref_colnames[1])])
          end_regime <- max(sub_aln_dt[(end_indexes[(i-1)]+1):end_indexes[i], get(start_end_ref_colnames[2])])
          
          start_regime_vec <- c(start_regime_vec, start_regime)
          end_regime_vec <- c(end_regime_vec, end_regime)
          
        }
        
        out_rle_dt <- data.table::data.table(start_regime=start_regime_vec, end_regime=end_regime_vec)
        out_rle_dt[, tmp_ratio_name:=rle_combination$values]
        colnames(out_rle_dt)[3] <- ratio_name
        out_rle_dt[, length:=end_regime-start_regime]
        
        by_combo_representation <- foreach::foreach(combo = unique(out_rle_dt[[ratio_name]]))%do%{
          which_row <- out_rle_dt[[ratio_name]]==combo
          total_this_combo <- sum(out_rle_dt$length[which_row])
          return(total_this_combo)
        }
        count_dt <- data.table::data.table(unique(out_rle_dt[[ratio_name]]), bp_length=unlist(by_combo_representation))
        count_dt$percentage <- count_dt$bp_length/sum(count_dt$bp_length)*100
        colnames(count_dt)[1] <- ratio_name
        
        return(list(count_dt=count_dt, run_length_encoding=out_rle_dt))
      }
      doParallel::stopImplicitCluster()
      
      names(per_chromo_list) <- chromo
      
      full_count_temp_dt <- data.table::rbindlist(lapply(per_chromo_list, function(chr_list){
        return(chr_list$count_dt)
      }))
      
      unique_ratios <- unique(full_count_temp_dt[[ratio_name]])
      total_count_list <- foreach::foreach(combo = unique_ratios)%do%{
        which_row <- full_count_temp_dt[[ratio_name]]==combo
        return(sum(full_count_temp_dt$bp_length[which_row]))
      }
      
      total_prog_count_dt <- data.table::data.table(bp_length=unlist(total_count_list))
      total_prog_count_dt$percentage <- total_prog_count_dt$bp_length/sum(total_prog_count_dt$bp_length)*100
      total_prog_count_dt[[ratio_name]] <- unique_ratios
      total_prog_count_dt <- total_prog_count_dt[order(total_prog_count_dt[[ratio_name]]),]
      
      total_name <- paste0("total_", prog)
      per_chromo_list[[total_name]] <- total_prog_count_dt
      
      return(per_chromo_list)
    }
    
    names(summary_per_prog_list) <- progenitors
    
    list_total_dts <- lapply(summary_per_prog_list, function(prog){
      total_name <- grep("total_", names(prog), value=TRUE)
      return(prog[[total_name]])
    })
    
    tmp_sample_out_dt <- data.table::rbindlist(list_total_dts)
    sum_list <- foreach::foreach(combo=unique(tmp_sample_out_dt[[ratio_name]]))%do%{
      
      which_row <- tmp_sample_out_dt[[ratio_name]]==combo
      return(sum(tmp_sample_out_dt$bp_length[which_row]))
      
    }
    smp_total_dt <- data.table::data.table(count=unlist(sum_list))
    smp_total_dt$percentage <- smp_total_dt$count/sum(smp_total_dt$count)*100
    smp_total_dt[[ratio_name]] <- unique(tmp_sample_out_dt[[ratio_name]])
    smp_total_dt <- smp_total_dt[order(smp_total_dt[[ratio_name]]),]
    
    summary_per_prog_list[["total_summary_dt"]] <- smp_total_dt
    
    return(summary_per_prog_list)
  }
  
  names(alignment_summary_list) <- polyploid_samples
  return(alignment_summary_list)
  
}


utils::globalVariables(c("combination", "combo", "tmp_ratio_name"))