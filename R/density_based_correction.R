#' Get densities from concordant bins.
#'
#' Get the two-dimensional (counts and GC content) kernel density estimation per copy number using bins involved in a uniquely overlapping and concordant anchor
#'
#' @param alignment A heal alignment object created with get_heal_alignment().
#' @param heal_list Output list in heal format (such as output from count_heal_data()).
#' @param n_threads Number of threads to use ('1' by default).
#' @param prog_ploidy Ploidy of the progenitors (Assumed to be equal. '2' by default).
#' @param n_points The n parameter for MASS::kde2d(). "Number of grid points in each direction. Can be scalar or a length-2 integer vector". Default is 1000.
#' @param normalize Use normalize count values (either FALSE, 'local' or 'manual'. FALSE by default). The 'local' and 'manual' entries have the same meaning as for 'method' in 'get_copy_number()'.  
#' @param average Average measure used normalize each bin count ('median' or 'mean'. 'median' by default).
#' @param average_list A list with one element per sample with each containing one element per subgenome with values used for normalization (same as in 'get_copy_number()').
#'
#' @return A list with one element per polyploid sample. Each sample has a list with one density distribution for each copy number from 0 to the ploidy of the progenitors times the number of progenitors.
#' @export
#'
get_concordant_density <- function(alignment, heal_list, n_threads = 1, prog_ploidy = 2, n_points = 1000, normalize = FALSE, average = "median", average_list = FALSE) {
  
  polyploid_samples <- names(alignment)
  progenitors <- names(heal_list)
  total_ploidy <- length(progenitors) * prog_ploidy
  
  if(normalize == "local"){
    local_averages <- get_sample_stats(heal_list, method = paste0("local_", average))
    
  }else if (normalize == "manual" & average_list != FALSE){
    if (!identical(names(average_list), polyploid_samples)){
      stop("Invalid average_list input. When method == 'manual' you must provide an average_list with values of median for each sample subgenome.")
    }
  }else if (normalize == "manual" & average_list == FALSE){
    stop("No 'average_list' provided. When method == 'manual' you must provide an 'average_list' with values of median for each sample subgenome.")
  }
  
  
  doParallel::registerDoParallel(n_threads)
  cn_count_type_dt_list <- foreach::foreach(smp = polyploid_samples) %dopar% {
    
    concordant_and_unique <- alignment[[smp]]$status == "concordant" & alignment[[smp]]$method == "unique"

    dt_per_prog_list <- foreach::foreach(prog = progenitors) %do% {
      
      # Names of columns of interest
      bin_index_col <- paste0("bin_index_", prog)
      
      # Extract index of bins which are concordant and unique
      cc_and_u_bins <- apply(alignment[[smp]][concordant_and_unique, ..bin_index_col], 1, function(row) {
        strsplit(row, ",")
      })

      cc_and_u_rows <- as.numeric(unique(unlist(cc_and_u_bins)))
      
      # Make a data table of only cc and u bins and their count values using "merge"
      chr_vec <- heal_list[[prog]]$CN$chr[cc_and_u_rows]
      start_vec <- heal_list[[prog]]$CN$start[cc_and_u_rows]
      cc_and_u_dt <- data.table::data.table(chr = chr_vec, start = start_vec)
      data.table::setkey(cc_and_u_dt, chr, start)
      cc_and_u_dt <- merge(heal_list[[prog]]$bins, cc_and_u_dt, by = c("chr", "start"))
      col_keep <- c("chr", "start", "gc_content", smp)
      cc_and_u_dt <- cc_and_u_dt[, ..col_keep]
      data.table::setkey(cc_and_u_dt, chr, start)
      colnames(cc_and_u_dt) <- c("chr", "start", "gc_content", "count") # rename to avoid overlapp with CN
      cc_and_u_dt <- merge(heal_list[[prog]]$CN, cc_and_u_dt, by = c("chr", "start", "gc_content"))
      col_keep <- c("chr", "start", "gc_content", "count", smp)
      cc_and_u_dt <- cc_and_u_dt[, ..col_keep]
      colnames(cc_and_u_dt) <- c("chr", "start", "gc_content", "count", "cn")

      cc_and_u_dt <- cc_and_u_dt[!is.na(cc_and_u_dt$count), ]
      
      if(normalize == "local"){
        cc_and_u_dt$count <- cc_and_u_dt$count / local_averages[[smp]][[prog]] * prog_ploidy
        
      }else if(normalize == "manual"){
        cc_and_u_dt$count <- cc_and_u_dt$count / average_list[[smp]][[prog]] * prog_ploidy
        
      }
      
      return(cc_and_u_dt)
    }

    total_dt <- data.table::rbindlist(dt_per_prog_list)
    key_colname <- c("chr", "start", "gc_content", "cn")
    data.table::setkeyv(total_dt, key_colname)

    return(total_dt)
  }
  doParallel::stopImplicitCluster()

  names(cn_count_type_dt_list) <- polyploid_samples

  doParallel::registerDoParallel(n_threads)
  density_list <- foreach::foreach(smp = polyploid_samples) %do% {
    cn_groups <- as.vector(na.omit(unique(cn_count_type_dt_list[[smp]]$cn)))
    
    density_by_cn_list <- foreach::foreach(cn = cn_groups) %do% {
      which_rows <- cn_count_type_dt_list[[smp]]$cn == cn & !is.na(cn_count_type_dt_list[[smp]]$cn)

      density <- MASS::kde2d(cn_count_type_dt_list[[smp]]$gc_content[which_rows], cn_count_type_dt_list[[smp]]$count[which_rows], n = n_points)

      return(density)
    }
    names(density_by_cn_list) <- paste0("cn_", cn_groups)
    return(density_by_cn_list)
  }
  doParallel::stopImplicitCluster()
  names(density_list) <- polyploid_samples

  return(density_list)
}

#' Re-evaluate the copy number at discordant anchors
#' Goes over all discordant anchor sets and verifies if some of their overlapping bins can be changed based on density. If a change leads to a concordant combination of copy number for the set of progenitors then it is accepted.
#'
#' @param densities A list of densities for each polyploid sample (output from get_concordant_density()).
#' @param alignment A heal alignment object created with get_heal_alignment().
#' @param heal_list Output list in heal format (such as output from count_heal_data()).
#' @param n_threads Number of threads to use ('1' by default).
#' @param prog_ploidy Ploidy of the progenitors (Assumed to be equal. '2' by default).
#'
#' @return A heal alignment object with the re-evaluated copy number based on empirical densities. The status column indicates what correction procedure (if any) was applied
#' @export
#'
#' @importFrom data.table :=
correct_cn_with_density <- function(densities, alignment, heal_list, n_threads, prog_ploidy = 2) {
  
  polyploid_samples <- names(densities)
  progenitors <- names(heal_list)

  corrected_alignment <- foreach::foreach(smp = polyploid_samples) %do% {
    
    corrected_aln_dt <- alignment[[smp]]
    which_discord <- which(corrected_aln_dt$status == "discordant")

    doParallel::registerDoParallel(n_threads)
    corrected_discordant_list <- foreach::foreach(a = which_discord) %dopar% {
      
      row <- corrected_aln_dt[a, ]
      
      valid_cn_at_anchor_list <- foreach::foreach(prog = progenitors) %do% {
        
        cn_col_name <- paste0("cn_", prog)
        bin_col_name <- paste0("bin_index_", prog)
        bindex_vec <- as.numeric(unlist(strsplit(row[[bin_col_name]], ",")))

        chr_vec <- heal_list[[prog]]$CN$chr[bindex_vec]
        start_vec <- heal_list[[prog]]$CN$start[bindex_vec]

        merged_dt <- merge(heal_list[[prog]]$bins, data.table::data.table(chr = chr_vec, start = start_vec), by = c("chr", "start"))
        col_keep <- c("gc_content", smp)
        merged_dt <- merged_dt[, ..col_keep]
        merged_dt$cn <- rep(row[[cn_col_name]], length(bindex_vec))

        valid_cn_list <- foreach::foreach(i = 1:nrow(merged_dt)) %do% {
          density_at_cn_list <- foreach::foreach(cn = names(densities[[smp]])) %do% {
            x_idx <- which.min(abs(densities[[smp]][[cn]]$x - merged_dt$gc_content[i]))
            y_idx <- which.min(abs(densities[[smp]][[cn]]$y - merged_dt[[smp]][i]))
            dens_bin <- densities[[smp]][[cn]]$z[x_idx, y_idx]
            return(dens_bin)
          }
          names(density_at_cn_list) <- names(densities[[smp]])
          density_at_cn <- unlist(density_at_cn_list)

          current_cn_name <- paste0("cn_", merged_dt$cn[i])

          if (length(intersect(current_cn_name, names(density_at_cn_list))) == 0) {
            which_best <- which.max(density_at_cn_list)
            valid_cn <- as.numeric(gsub("cn_", "", names(density_at_cn_list)[which_best]))
            densities_valid <- unlist(density_at_cn_list[which_best])
            return(data.table::data.table(cn = valid_cn, densities = densities_valid))
          } else {
            which_better_or_equal <- density_at_cn >= density_at_cn[current_cn_name]

            valid_cn <- as.numeric(gsub("cn_", "", names(density_at_cn_list)[which_better_or_equal]))
            densities_valid <- unlist(density_at_cn_list[which_better_or_equal])
            return(data.table::data.table(cn = valid_cn, densities = densities_valid))
          }
        }

        possible_cn <- unique(unlist(lapply(valid_cn_list, function(dt) dt$cn)))

        best_density <- foreach::foreach(cn = possible_cn) %do% {
          return(max(unlist(lapply(valid_cn_list, function(dt) {
            return(dt$densities[dt$cn == cn])
          }))))
        }

        return(data.table::data.table(cn = possible_cn, max_density = unlist(best_density)))
      }

      names(valid_cn_at_anchor_list) <- progenitors

      cn_only_list <- lapply(valid_cn_at_anchor_list, function(dt) dt$cn)
      all_combinations <- expand.grid(cn_only_list)
      which_concordant <- rowSums(all_combinations) == prog_ploidy * length(progenitors)

      concordant_combinations <- all_combinations[which_concordant, ]

      if (sum(which_concordant) > 1) {
        sum_density_list <- apply(concordant_combinations, 1, function(row) {
          density_sum <- 0
          for (i in 1:length(row)) {
            prog <- names(row)[i]
            which_cn_row <- valid_cn_at_anchor_list[[prog]]$cn == row[i]
            density_sum <- density_sum + valid_cn_at_anchor_list[[prog]]$max_density[which_cn_row]
          }

          return(density_sum)
        })

        which_max_density_combination <- which.max(unlist(sum_density_list))
        replacement_dt <- concordant_combinations[which_max_density_combination, ]
        colnames(replacement_dt) <- paste0("cn_", colnames(replacement_dt))

        for (col_name in colnames(replacement_dt)) {
          data.table::set(corrected_aln_dt, i = a, j = col_name, value = replacement_dt[[col_name]])
        }
        corrected_aln_dt[a, status := "corrected_multiple"]
      } else if (sum(which_concordant) == 1) {
        replacement_dt <- concordant_combinations
        colnames(replacement_dt) <- paste0("cn_", colnames(replacement_dt))

        for (col_name in colnames(replacement_dt)) {
          data.table::set(corrected_aln_dt, i = a, j = col_name, value = replacement_dt[[col_name]])
        }
        corrected_aln_dt[a, status := "corrected_single"]
      }
      return(corrected_aln_dt[a, ])
    }
    doParallel::stopImplicitCluster()

    merge_dt <- data.table::rbindlist(corrected_discordant_list)
    merge_dt <- merge(corrected_aln_dt, merge_dt, by = paste0("id_", progenitors), all.x = TRUE, suffixes = c("_original", "_corrected"))

    col_names_change <- c(paste0("cn_", progenitors), "status")
    for (col_name in col_names_change) {
      original_colname <- paste0(col_name, "_original")
      corrected_colname <- paste0(col_name, "_corrected")
      merge_dt[!is.na(get(corrected_colname)), (original_colname) := get(corrected_colname)]
    }

    colnames_remove <- grep("_corrected$", colnames(merge_dt), value = TRUE)
    merge_dt[, (colnames_remove) := NULL]

    colnames_rename <- grep("_original$", colnames(merge_dt), value = TRUE)
    replacement <- sub("_original", "", colnames_rename)
    data.table::setnames(merge_dt, colnames_rename, replacement)

    return(merge_dt)
  }

  names(corrected_alignment) <- names(alignment)
  return(corrected_alignment)
}


utils::globalVariables(c("a", "col_keep", "cn", "i", "status", "bin_index_col", "which_relevant", "cn_col", "bin_start", "anchor_start", "anchor_end", "gene_id", "ref", "copy", "subgenome", "points", "progenitor", "counts", "is_allowed", "x", "y"))
