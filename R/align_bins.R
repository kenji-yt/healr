#' Get the genome pairs of the synHits file.
#'
#' @param syn_hits_list A list with the GENESPACE output files "synHits.txt.gz".
#' @param show_which_non_self Output a logical value being TRUE if the pair is not a self hit. Defaults is "TRUE".
#'
#' @return A list with one element per genome pair (including self). If show_which_non_self=TRUE the item is a logical value being TRUE if the pair is self hit, else it returns the base name of the file without suffix.
#'
ref_alt_pair <- function(syn_hits_list, show_which_non_self = TRUE) {
  out_list <- lapply(syn_hits_list, function(path) {
    filename <- basename(path)
    filename <- gsub(".synHits.txt.gz", "", filename)
    which_genomes <- unique(unlist(strsplit(filename, "_vs_")))
    if (show_which_non_self == TRUE) {
      return(length(which_genomes) != 1)
    } else {
      paste(filename)
    }
  })
  return(unlist(out_list))
}

#' Parse genespace results directory to extract syntenic hits information.
#' Gets a list with one data table per pair of progenitors. Each data table contains all syntenic anchor pairs (syntenicHits.txt files from GENESPACE).
#'
#' @param genespace_dir Path to a directory containing the syntenicHits output directory for GENESPACE.
#'
#' @return A list with one data table per pair of progenitors. Each data table contains all syntenic anchor pairs (syntenicHits.txt files from GENESPACE).
#'
parse_genespace_input <- function(genespace_dir) {
  # Syntenic hits
  all_syn_lists <- list.files(genespace_dir, pattern = "synHits.txt.gz$", recursive = TRUE, full.names = TRUE)
  which_keep <- ref_alt_pair(all_syn_lists, show_which_non_self = TRUE)
  syn_hits_list <- list()
  for (i in 1:length(all_syn_lists)) {
    if (which_keep[i] == TRUE) syn_hits_list <- append(syn_hits_list, all_syn_lists[[i]])
  }
  syn_hits_dt_list <- lapply(syn_hits_list, function(synHits_txt) {
    return(data.table::as.data.table(utils::read.table(gzfile(synHits_txt[[1]]), header = TRUE, sep = "\t")))
  })

  names(syn_hits_dt_list) <- ref_alt_pair(syn_hits_list, show_which_non_self = FALSE)


  return(syn_hits = syn_hits_dt_list)
}


#' get_conserved_hits
#'
#' @param genespace_dir Path to a directory containing the syntenicHits output directory for GENESPACE.
#'
#' @return A data table with one row for each anchor set present in all progenitors.
#'
#' @importFrom data.table .SD
get_conserved_anchors <- function(genespace_dir) {
  syn_hits_list <- parse_genespace_input(genespace_dir)

  # Pick on genome as reference and keep only data tables containing this genome.
  genomes_per_hits_dt_list <- strsplit(names(syn_hits_list), "_vs_")
  reference_genome <- unlist(genomes_per_hits_dt_list)[1]
  which_has_ref <- which(sapply(genomes_per_hits_dt_list, function(x) reference_genome %in% x))

  genomes <- c(syn_hits_list[which_has_ref][[1]]$genome1[1], syn_hits_list[which_has_ref][[1]]$genome2[1])
  ref_index <- which(genomes == reference_genome)
  ref_chr <- paste0("chr", ref_index)
  ref_start <- paste0("start", ref_index)
  ref_end <- paste0("end", ref_index)
  ref_id <- paste0("id", ref_index)

  merge_dt <- syn_hits_list[which_has_ref][[1]][, .SD, .SDcols = c(ref_chr, ref_start, ref_end, ref_id)]
  merge_dt <- merge_dt[syn_hits_list[which_has_ref][[1]]$isAnchor == TRUE, ]
  colnames(merge_dt) <- c(paste0("chr_", reference_genome), paste0("start_", reference_genome), paste0("end_", reference_genome), paste0("id_", reference_genome))

  for (hits_dt in syn_hits_list[which_has_ref]) {
    anchors_dt <- hits_dt[hits_dt$isAnchor == TRUE, ]
    genome_1 <- anchors_dt$genome1[1]
    genome_2 <- anchors_dt$genome2[1]

    colnames(anchors_dt)[colnames(anchors_dt) == c("chr1")] <- paste0("chr_", genome_1)
    colnames(anchors_dt)[colnames(anchors_dt) == c("start1")] <- paste0("start_", genome_1)
    colnames(anchors_dt)[colnames(anchors_dt) == c("end1")] <- paste0("end_", genome_1)
    colnames(anchors_dt)[colnames(anchors_dt) == c("id1")] <- paste0("id_", genome_1)

    colnames(anchors_dt)[colnames(anchors_dt) == c("chr2")] <- paste0("chr_", genome_2)
    colnames(anchors_dt)[colnames(anchors_dt) == c("start2")] <- paste0("start_", genome_2)
    colnames(anchors_dt)[colnames(anchors_dt) == c("end2")] <- paste0("end_", genome_2)
    colnames(anchors_dt)[colnames(anchors_dt) == c("id2")] <- paste0("id_", genome_2)

    merge_dt <- merge(merge_dt, anchors_dt, by = intersect(names(merge_dt), names(anchors_dt)))

    which_relevant <- which(colnames(merge_dt) %in% c(
      paste0(c("chr_", "start_", "end_", "id_"), genome_1),
      paste0(c("chr_", "start_", "end_", "id_"), genome_2)
    ))

    merge_dt <- merge_dt[, ..which_relevant]
  }

  return(merge_dt)
}



#' Get the bins overlapping with each anchor gene.
#'
#' @param anchors_dt A data table with one row for each anchor set present in all progenitors (such as output from get_conserved_anchors()).
#' @param heal_list Output list in heal format (such as output from count_heal_data()).
#' @param genespace_dir Path to a directory containing the syntenicHits output directory for GENESPACE.
#' @param n_threads Number of threads to use ('1' by default).
#'
#' @return A list containing one data table per progenitor with anchors and overlapping bins.
#'
map_bins_to_anchors <- function(anchors_dt, heal_list, genespace_dir, n_threads=1) {
  syn_hits_list <- parse_genespace_input(genespace_dir)

  sample_names <- unlist(lapply(heal_list, function(prog) {
    setdiff(colnames(prog$bins), c("chr", "start", "end", "mappability", "gc_content"))
  }))
  polyploid_samples <- names(table(sample_names))[table(sample_names) == 2]

  anchors_dt <- get_conserved_anchors(genespace_dir)

  progenitors <- names(heal_list)

  syn_anchors_to_bin_dt_list <- list()
  for (prog in progenitors) {
    chr <- paste0("chr_", prog)
    start <- paste0("start_", prog)
    end <- paste0("end_", prog)
    id <- paste0("id_", prog)

    gnm_anchors_ranges <- GenomicRanges::GRanges(
      seqnames = anchors_dt[[chr]],
      ranges = IRanges::IRanges(
        start = anchors_dt[[start]],
        end = anchors_dt[[end]]
      ),
      gene_id = anchors_dt[[id]]
    )

    gnm_bin_ranges <- GenomicRanges::GRanges(
      seqnames = heal_list[[prog]]$CN$chr,
      ranges = IRanges::IRanges(
        start = heal_list[[prog]]$CN$start,
        end = heal_list[[prog]]$CN$end
      )
    )

    gnm_overlap_bins_anchors <- GenomicRanges::findOverlaps(query = gnm_bin_ranges, subject = gnm_anchors_ranges)

    dt_overlap_bins_anchors <- data.table::as.data.table(gnm_overlap_bins_anchors) # Because queryHits does not work..?

    gnm_overlap_widths <- IRanges::width(IRanges::pintersect(gnm_bin_ranges[dt_overlap_bins_anchors$queryHits], gnm_anchors_ranges[dt_overlap_bins_anchors$subjectHits]))

    gnm_overlap_dt <- data.table::data.table(
      bin_index = dt_overlap_bins_anchors$queryHits,
      anchor_index = dt_overlap_bins_anchors$subjectHits,
      bin_start = IRanges::start(gnm_bin_ranges[dt_overlap_bins_anchors$queryHits]),
      anchor_start = IRanges::start(gnm_anchors_ranges[dt_overlap_bins_anchors$subjectHits]),
      anchor_end = IRanges::end(gnm_anchors_ranges[dt_overlap_bins_anchors$subjectHits]),
      gene_id = gnm_anchors_ranges$gene_id[dt_overlap_bins_anchors$subjectHits],
      gc_content = heal_list[[prog]]$CN$gc_content[dt_overlap_bins_anchors$queryHits],
      overlap_width = gnm_overlap_widths
    )

    doParallel::registerDoParallel(n_threads)
    cn_per_sample_list <- foreach::foreach(smp = polyploid_samples) %dopar% {
      return(heal_list[[prog]]$CN[[smp]][dt_overlap_bins_anchors$queryHits])
    }
    doParallel::stopImplicitCluster()

    cn_dt <- data.table::as.data.table(cn_per_sample_list)
    colnames(cn_dt) <- polyploid_samples

    cn_anchors_and_bins_dt <- cbind(gnm_overlap_dt, cn_dt)
    
    data.table::setkey(cn_anchors_and_bins_dt, bin_start, anchor_start, anchor_end, gene_id, gc_content)
    
    cn_anchors_and_bins_dt <- merge(anchors_dt, cn_anchors_and_bins_dt, by.x = paste0("id_", prog), by.y = "gene_id")
    
    syn_anchors_to_bin_dt_list[[prog]] <- cn_anchors_and_bins_dt
  }

  return(syn_anchors_to_bin_dt_list)
}


#' Create heal alignment object.
#'
#' @param heal_list Output list in heal format (such as output from count_heal_data()).
#' @param genespace_dir Path to a directory containing the syntenicHits output directory for GENESPACE.
#' @param n_threads Number of threads to use ('1' by default).
#' @param prog_ploidy Ploidy of the progenitors (Assumed to be equal. '2' by default).
#'
#' @return A list with one data table per sample. Each data table has one row for each anchor set present in all progenitors. The bins overlapping each anchor are given along with a consensus copy number at the anchor based on these bins.
#' @export
#' #'
get_heal_alignment <- function(heal_list, genespace_dir, n_threads = 1, prog_ploidy = 2) {

  cn_exist <- sum(names(heal_list[[1]]) == "CN") != 0
  if (cn_exist != TRUE) {
    stop("No CN data. Exiting...")
  }

  sample_names <- unlist(lapply(heal_list, function(prog) {
    setdiff(colnames(prog$bins), c("chr", "start", "end", "mappability", "gc_content"))
  }))
  polyploid_samples <- names(table(sample_names))[table(sample_names) == 2]
  progenitors <- names(heal_list)

  total_ploidy <- length(progenitors) * prog_ploidy

  # Get data table of anchors present in all subgenomes
  anchors_dt <- get_conserved_anchors(genespace_dir)
  # Get a list with one data table per subgenome containing anchors and their overlapping bins
  cn_anchors_and_bins <- map_bins_to_anchors(anchors_dt, heal_list, genespace_dir, n_threads)
#
#   col_to_merge <- colnames(cn_anchors_and_bins[[progenitors[1]]])
#   col_to_merge <- col_to_merge[grepl(paste0("^", progenitors[1], "_"), col_to_merge )]
#
#   merged_dt <- Reduce(function(x, y) merge(x, y, by=col_to_merge, all = TRUE), cn_anchors_and_bins)

  # Solution is to use join to merge the anchors dt with the cn_anchors_and_bins data tables.
  # Once merged (can probably merge all subgenome dt into one) we can go through that row by row and assign a CN.
  # Sweet!
  cn_alignment_list <- foreach::foreach(smp = polyploid_samples) %do% {

    assigned_and_unassigned_list <- foreach::foreach(prog = progenitors) %do% {
      
      prog_id <- paste0("id_", prog)
      # We can extract anchors which overlap multiple bins and process these individually
      multiple_overlap_anchors <- table(cn_anchors_and_bins[[prog]][[prog_id]])[table(cn_anchors_and_bins[[prog]][[prog_id]])>1]

      # First make a dt for the uniquely overlapping anchors
      unique_assigned_anchors <- setdiff(cn_anchors_and_bins[[prog]][[prog_id]], names(multiple_overlap_anchors))
      
      assigned_anchors_dt <- cn_anchors_and_bins[[prog]][get(prog_id) %in% unique_assigned_anchors]
      
      # remove superfluous columns
      cols_to_remove <- c("anchor_index", "bin_start", "anchor_start", "anchor_end", "gc_content", "overlap_width")
      assigned_anchors_dt[, (cols_to_remove) := NULL]
      
      # add necessary columns
      colnames(assigned_anchors_dt)[colnames(assigned_anchors_dt)==smp] <- paste0("cn_", prog)
      colnames(assigned_anchors_dt)[colnames(assigned_anchors_dt)=="bin_index"] <- paste0("bin_index_", prog)
      assigned_anchors_dt$method <- rep("unique", nrow(assigned_anchors_dt))
      
      # Now we can see which multiple overlapping ones have a unique 
      doParallel::registerDoParallel(n_threads)
      anchor_method_list <- foreach::foreach(i = 1:length(multiple_overlap_anchors)) %dopar% {
        
        anchor_id <- names(multiple_overlap_anchors[i])
        sub_dt <- cn_anchors_and_bins[[prog]][cn_anchors_and_bins[[prog]][[prog_id]]==anchor_id]
        cn_vec <- sub_dt[[smp]]
        
        if(length(unique(cn_vec))>1){
          method <- "multiple"
          overlapp_vec <- sub_dt$overlap_width
          bin_index <- sub_dt$bin_index
          method_cn_list <- list(method = method, cn = unique(cn_vec), overlapp = overlapp_vec, anchor = anchor_id, bin_index = bin_index)
        }else{
          method <- "unique"
          bin_index <- sub_dt$bin_index
          method_cn_list <- list(method=method, cn=unique(cn_vec), anchor = anchor_id, bin_index = bin_index)
        }
        
        return(method_cn_list)
      }
      doParallel::stopImplicitCluster()
      
      # We find them and add them to the assigned_anchors_dt
      unique_prog_list <- foreach::foreach(i = 1:length(anchor_method_list))%dopar%{
        
        if(anchor_method_list[[i]]$method=="unique"){
          
          dt <- data.table::data.table(id = anchor_method_list[[i]]$anchor, cn = anchor_method_list[[i]]$cn, method = "unique", bin_index = anchor_method_list[[i]]$bin_index)
          colnames(dt)[colnames(dt)=="id"] <- prog_id
          colnames(dt)[colnames(dt)=="cn"] <- paste0("cn_", prog)
          colnames(dt)[colnames(dt)=="bin_index"] <- paste0("bin_index_", prog)
          
          return(dt)
        }
      }
      doParallel::stopImplicitCluster()
      
      unique_prog_dt <- data.table::rbindlist(unique_prog_list)
      assigned_add_on_dt<- merge(unique_prog_dt, anchors_dt, by = prog_id)
      
      assigned_anchors_dt <- data.table::rbindlist(list(assigned_anchors_dt, assigned_add_on_dt), use.names = TRUE)
      
      # Now we format the non-unique ones and return then for later assignment.
      multiple_prog_list <- Filter(function(x) x$method == "multiple", anchor_method_list)
      
      return(list(assigned_dt = assigned_anchors_dt, unassigned_list = multiple_prog_list))
    }
    
    doParallel::registerDoParallel(n_threads)
    cn_per_anchor_pair_list <- foreach::foreach(i = 1:nrow(merged_dt)) %dopar% {

      # cn_at_anchor_dt_list <- foreach::foreach(prog = progenitors) %do% {
      #
      #   anchor_gene <- anchors_dt[[paste0("id_", prog)]][i]
      #
      #   which_rows <- cn_anchors_and_bins[[prog]][[paste0("id_", prog)]] == anchor_gene
      #
      #   cn_at_anchor_dt <- cn_anchors_and_bins[[prog]][which_rows, ]
      #
      #   return(cn_at_anchor_dt)
      # }
      # names(cn_at_anchor_dt_list) <- progenitors

      # # If the anchor overlaps no bin (eg. bin filtered out)
      # if (sum(lapply(cn_at_anchor_dt_list, nrow) == 0) != 0) {
      #   return()
      # } else {
      unique_cn <- foreach::foreach(prog = progenitors) {
        merged_dt[i,get(paste0(""))]
        unique(dt[[smp]])
      })

        bin_indexes <- lapply(cn_at_anchor_dt_list, function(dt) {
          paste(dt$bin_index, collapse = ",")
        })

        bin_index_dt <- data.table::as.data.table(bin_indexes)
        colnames(bin_index_dt) <- paste0("bin_index_", colnames(bin_index_dt))

        # Unique CN in all genomes
        if (sum(unlist(lapply(unique_cn, length))) == length(unique_cn)) {
          unique_cn_dt <- data.table::as.data.table(unique_cn)
          colnames(unique_cn_dt) <- paste0("cn_", colnames(unique_cn_dt))

          output_dt <- cbind(anchors_dt[i, ], unique_cn_dt, bin_index_dt, data.table::data.table(method = "unique"))
          return(output_dt)

          # Find most concordant CN combination (heuristic solution)
        } else {
          cartesian_product <- data.table::as.data.table(expand.grid(unique_cn))
          most_concordant <- which(abs(rowSums(cartesian_product) - total_ploidy) == min(abs(rowSums(cartesian_product) - total_ploidy)))

          concordant_combinations <- cartesian_product[most_concordant, ]

          # if multiple equally concordant combinations, pick the one with most overlap
          # Note that we consider copy number of same value at start and end of the gene as equivalent
          if (nrow(concordant_combinations) > 1) {
            overlaps <- apply(concordant_combinations, 1, function(row) {
              overlaps <- foreach::foreach(prog = names(row)) %do% {
                which_row <- cn_at_anchor_dt_list[[prog]][[smp]] == row[[prog]]
                overlap <- sum(cn_at_anchor_dt_list[[prog]]$overlap_width[which_row])
                return(overlap)
              }
              sum(unlist(overlaps))
            })

            unique_cn_dt <- concordant_combinations[which.max(overlaps)]
            colnames(unique_cn_dt) <- paste0("cn_", colnames(unique_cn_dt))

            output_dt <- cbind(anchors_dt[i, ], unique_cn_dt, bin_index_dt, data.table::data.table(method = "multiple_overlap"))
            return(output_dt)
          } else {
            unique_cn_dt <- concordant_combinations
            colnames(unique_cn_dt) <- paste0("cn_", colnames(unique_cn_dt))

            output_dt <- cbind(anchors_dt[i, ], unique_cn_dt, bin_index_dt, data.table::data.table(method = "multiple_unique"))
            return(output_dt)
          }
        }
      }
    }
    doParallel::stopImplicitCluster()

    alignment <- data.table::rbindlist(cn_per_anchor_pair_list)

    alignment$status <- rep("concordant", nrow(alignment))

    cn_col <- grep("cn_", colnames(alignment))
    discordant_cols <- rowSums(alignment[, ..cn_col]) != prog_ploidy * length(heal_list)
    alignment$status[discordant_cols] <- "discordant"
    return(alignment)
  }

  names(cn_alignment_list) <- polyploid_samples
  return(cn_alignment_list)
}
