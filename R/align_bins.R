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
#' @param genespace_dir Path to a directory containing the syntenicHits output directory from GENESPACE.
#' @param n_threads Number of threads to use ('1' by default).
#' @param prog_ploidy Ploidy of the progenitors (Assumed to be equal. '2' by default).
#'
#' @return A list with one data table per sample. Each data table has one row for each anchor set present in all progenitors. The bins overlapping each anchor are given along with a consensus copy number at the anchor based on these bins.
#' @export
#'
get_heal_alignment <- function(heal_list, genespace_dir, n_threads = 1, prog_ploidy = 2) {

  cn_exist <- sum(names(heal_list[[1]]) == "CN") != 0
  if (cn_exist != TRUE) {
    stop("No CN data. Exiting...")
  }

  bin_size <- heal_list$A.arenosa$bins$end[1] - heal_list$A.arenosa$bins$start[1]
  
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

  # Make alignment for each sample
  cn_alignment_list <- foreach::foreach(smp = polyploid_samples) %do% {
    
    assigned_list <- foreach::foreach(prog = progenitors) %do% {
      
      prog_id_col <- paste0("id_", prog)
      # We can extract anchors which overlap multiple bins. We process these individually later.
      multiple_overlap_anchors <- table(cn_anchors_and_bins[[prog]][[prog_id_col]])[table(cn_anchors_and_bins[[prog]][[prog_id_col]])>1]

      # First make a dt for the uniquely overlapping anchors
      unique_assigned_anchors <- setdiff(cn_anchors_and_bins[[prog]][[prog_id_col]], names(multiple_overlap_anchors))
      
      assigned_anchors_dt <- cn_anchors_and_bins[[prog]][get(prog_id_col) %in% unique_assigned_anchors]
      
      # Get name of other samples to remove from assigned_anchors_dt
      non_focal_samples <- setdiff(polyploid_samples, smp)
      
      # remove superfluous columns
      cols_to_remove <- c("anchor_index", "bin_start", "anchor_start", "anchor_end", "gc_content", "overlap_width", non_focal_samples)
      assigned_anchors_dt[, (cols_to_remove) := NULL]
      # add desired columns
      colnames(assigned_anchors_dt)[colnames(assigned_anchors_dt)==smp] <- paste0("cn_", prog)
      colnames(assigned_anchors_dt)[colnames(assigned_anchors_dt)=="bin_index"] <- paste0("bin_index_", prog)
      assigned_anchors_dt[[paste0(prog, "_method")]] <- rep("unique", nrow(assigned_anchors_dt))
      
      
      
      # Assign CN to anchors overlapping multiple based on overlap length
      doParallel::registerDoParallel(n_threads)
      multiple_cn_anchors <- foreach::foreach(i = 1:length(multiple_overlap_anchors)) %do% {
        
        anchor_id <- names(multiple_overlap_anchors[i])
        sub_dt <- cn_anchors_and_bins[[prog]][cn_anchors_and_bins[[prog]][[prog_id_col]]==anchor_id]
        cn_vec <- sub_dt[[smp]]
        
        # if more than one unique overlapping CN
        if(length(unique(cn_vec))>1){
          
          # Take the CN with largest overlap
          overlap_by_CN <- tapply(sub_dt$overlap_width, cn_vec, sum)
          most_overlap_cn <- as.numeric(names(which.max(overlap_by_CN)))
          bin_indexes <- paste(sub_dt$bin_index, collapse = ",")
          out_list <- list(anchor_id, most_overlap_cn, "multiple", bin_indexes)
          names(out_list) <- c(prog_id_col, paste0("cn_", prog), paste0(prog, "_method"), paste0("bin_index_", prog))
          
          return(out_list)
          
        }else{
          
          bin_indexes <- paste(sub_dt$bin_index, collapse = ",")
          out_list <- list(anchor_id, unique(cn_vec), "unique", bin_indexes)
          names(out_list) <- c(prog_id_col, paste0("cn_", prog), paste0(prog, "_method"), paste0("bin_index_", prog))
          
          return(out_list)
        
        }
      }
      doParallel::stopImplicitCluster()
      
      assigned_multiple_dt <- data.table::rbindlist(multiple_cn_anchors)
      
      assigned_add_on_dt <- merge(assigned_multiple_dt, anchors_dt, by = prog_id_col)
      
      assigned_anchors_dt <- data.table::rbindlist(list(assigned_anchors_dt, assigned_add_on_dt), use.names = TRUE)

      return(assigned_anchors_dt)
    }
    
    names(assigned_list) <- progenitors
    
    common_cols <- Reduce(intersect, lapply(assigned_list, names))
    alignment_draft <- Reduce(function(...) merge(..., by = common_cols, all = FALSE), assigned_list)
    
    method_cols <- grep("_method$", names(alignment_draft), value = TRUE)
    # Create new column "method" by comparing values row-wise
    alignment_draft[, method := ifelse(apply(.SD, 1, function(x) length(unique(x)) == 1), "unique", "multiple"), .SDcols = method_cols]
    alignment_draft[, (method_cols) := NULL]
    
    cn_cols <- grep("^cn_", names(alignment_draft), value = TRUE)
    alignment_draft[, status := ifelse(apply(.SD, 1, function(x) sum(x) == total_ploidy), "concordant", "discordant"), .SDcols = cn_cols]
    
    # Review discordant and multiple to see if another assignment can make them concordant
    mult_disc_rows <- which(alignment_draft$status == "discordant" & alignment_draft$method == "multiple")
      
    doParallel::registerDoParallel(n_threads)
    corrected_cn_list <- foreach::foreach(i = mult_disc_rows)%dopar%{
   
      discordant_row <- alignment_draft[i, ]
      
      # Get the CN and overlap at that anchor for each progenitor subgenome 
      cn_list <- list()
      overlap_list <- list()
      for(p in 1:length(progenitors)){

        prog <- progenitors[p]
        
        prog_id_col <- paste0("id_", prog)
        anchor_id <- alignment_draft[i, get(prog_id_col)]
        sub_dt <- cn_anchors_and_bins[[prog]][get(prog_id_col) %in% anchor_id]
        
        cn_vec <- sub_dt[[smp]]
        cn_list[[p]] <- unique(cn_vec)
        
        overlap_by_CN <- tapply(sub_dt$overlap_width, cn_vec, sum)
        overlap_list[[p]] <- overlap_by_CN
      
      }
      names(cn_list) <- progenitors
      names(overlap_list) <- progenitors

      combinations <- expand.grid(cn_list)
      
      # Keep only combinations without NA
      combinations <- combinations[rowSums(is.na(combinations))==0, ]
        
      # single
      if(sum(rowSums(combinations) == total_ploidy)==1){
        
        concord_combo <- combinations[rowSums(combinations) == total_ploidy,]
        return(concord_combo)
        
      # multiple concordant alternatives
      }else if(sum(rowSums(combinations) == total_ploidy)>1){
        
        # go through each concordant combination and find the one maximizing overlap
        concord_dt <- combinations[rowSums(combinations) == total_ploidy,]
        overlap_vec <- c()
        for(r in 1:nrow(concord_dt)){
          
          overlap <- 0
          for(prog in progenitors){
            cn_prog <- concord_dt[[prog]]
            overlap <- overlap + overlap_list[[prog]][[as.character(cn_prog)]]
          }
          overlap_vec <- c(overlap_vec, overlap)
        }
        
        return(concord_dt[which.max(overlap_vec),])
        
      # No concordant combination
      }else{
        return(NULL)
      }
    }
    doParallel::stopImplicitCluster()
  
    # Check which rows be replaced 
    to_replace <- !sapply(corrected_cn_list, is.null)
    
    if(sum(to_replace)!=0){
      rows_to_replace <- mult_disc_rows[to_replace]
      replacement_values <- corrected_cn_list[to_replace]
      cols_to_replace <- paste0("cn_", progenitors)
      
      # replace in the draft alignment
      for(i in 1:length(rows_to_replace)){
        alignment_draft[rows_to_replace[i],  (cols_to_replace) := as.list(replacement_values[[i]])]
        alignment_draft[rows_to_replace[i], (c("method", "status")) := as.list(c("multiple_concordant", "concordant"))]
      }
    }
   # for some reason it doesn't evaluate the data table after this loop otherwise 
    return(alignment_draft)
    
  }
  
  names(cn_alignment_list) <- polyploid_samples
  return(cn_alignment_list)
  
}
