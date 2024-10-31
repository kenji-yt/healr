ref_alt_pair <- function(syn_hits_list,show_which_non_self=TRUE){
  out_list <- lapply(syn_hits_list,function(path){
    filename <- basename(path)
    filename <- gsub(".synHits.txt.gz","",filename)
    which_genomes <- unique(unlist(strsplit(filename,"_vs_")))
    if(show_which_non_self==TRUE){
      return(length(which_genomes)!=1)
    }else{
      paste(filename)
    }
  })
  return(unlist(out_list))
}

parse_genespace_input <- function(genespace_dir){

  # Syntenic hits
  all_syn_lists <- list.files(genespace_dir,pattern = "synHits.txt.gz$",recursive = TRUE,full.names = TRUE)
  which_keep <- ref_alt_pair(all_syn_lists,show_which_non_self = TRUE)
  syn_hits_list <- list()
  for(i in 1:length(all_syn_lists)){
    if(which_keep[i]==TRUE)  syn_hits_list <- append(syn_hits_list,all_syn_lists[[i]])
  }
  syn_hits_dt_list <- lapply(syn_hits_list,data.table::fread)
  names(syn_hits_dt_list) <- ref_alt_pair(syn_hits_list,show_which_non_self = FALSE)


  # Block coordinates
  block_coord_path <- list.files(genespace_dir,pattern = "syntenicBlock_coordinates.csv",recursive = TRUE,full.names = TRUE)
  blk_coord_dt <- data.table::fread(block_coord_path)


  return(list(syn_hits=syn_hits_dt_list,block_coord=blk_coord_dt))
}

#' get_bins_to_map
#'
#' @param heal_list
#' @param blk_coord_dt
#'
#' @return
#'
#' @importFrom dplyr %>%
get_bins_to_map <- function(heal_list, blk_coord_dt){

  genomes <- unique(c(blk_coord_dt$genome1))

  bins_to_map_list <- foreach::foreach(gnm=genomes)%do%{


    current_blk_coord <- blk_coord_dt[blk_coord_dt$genome1==gnm & !grepl("selfBlk", blk_coord_dt$blkID), ]

    blk_ranges <- GenomicRanges::GRanges(seqnames = current_blk_coord$chr1,
                                        ranges = IRanges::IRanges(start = current_blk_coord$startBp1,
                                                                  end = current_blk_coord$endBp1),
                                        gene_id = current_blk_coord$blkID)

    bin_ranges <- GenomicRanges::GRanges(seqnames = heal_list[[gnm]]$CN$chr,
                                   ranges = IRanges::IRanges(start = heal_list[[gnm]]$CN$start,
                                                    end = heal_list[[gnm]]$CN$end))

    overlap_bins_blk <- GenomicRanges::findOverlaps(query = bin_ranges, subject = blk_ranges)

    dt_overlap_bins_blk <- data.table::as.data.table(overlap_bins_blk) #Because queryHits does not work..?

    overlap_widths <- IRanges::width(IRanges::pintersect(bin_ranges[dt_overlap_bins_blk$queryHits], blk_ranges[dt_overlap_bins_blk$subjectHits]))

    overlap_dt <- data.table::data.table(
      bin_index = dt_overlap_bins_blk$queryHits,
      blk_index = dt_overlap_bins_blk$subjectHits,
      bin_start = IRanges::start(bin_ranges[dt_overlap_bins_blk$queryHits]),
      bin_end = IRanges::end(bin_ranges[dt_overlap_bins_blk$queryHits]),
      bin_chr = heal_list[[gnm]]$CN$chr[dt_overlap_bins_blk$queryHits],
      blk_id = blk_ranges$gene_id[dt_overlap_bins_blk$subjectHits],
      overlap_width = overlap_widths
    )
#
#     overlap_dt <- overlap_dt %>%
#       dplyr::group_by(bin_index) %>%
#       dplyr::slice_max(order_by = overlap_width, with_ties = FALSE) %>%
#       dplyr::ungroup()

    return(overlap_dt)
  }
  names(bins_to_map_list) <- genomes
  return(bins_to_map_list)
}


#' get_conserved_hits
#'
#' @param genespace_dir
#'
#' @return
#'
get_conserved_anchors <- function(genespace_dir){

  map_dt_list <- parse_genespace_input(genespace_dir)
  syn_hits_list <- map_dt_list$syn_hits

  genomes_per_hits_dt_list <- strsplit(names(syn_hits_list),"_vs_")
  reference_genome <- unlist(genomes_per_hits_dt_list)[1]
  which_has_ref <- which(sapply(genomes_per_hits_dt_list, function(x) reference_genome %in% x))

  genomes <- c(syn_hits_list[which_has_ref][[1]]$genome1[1], syn_hits_list[which_has_ref][[1]]$genome2[1])
  ref_index <- which(genomes==reference_genome)
  ref_chr <- paste0("chr",ref_index)
  ref_start <- paste0("start",ref_index)
  ref_end <- paste0("end",ref_index)
  ref_id <- paste0("id",ref_index)

  merge_dt <- syn_hits_list[which_has_ref][[1]][,.SD, .SDcols = c(ref_chr, ref_start, ref_end, ref_id)]
  merge_dt <- merge_dt[syn_hits_list[which_has_ref][[1]]$isAnchor==TRUE, ]
  colnames(merge_dt) <- c(paste0("chr_",reference_genome), paste0("start_",reference_genome), paste0("end_",reference_genome), paste0("id_",reference_genome))

  for(hits_dt in syn_hits_list[which_has_ref]){

    anchors_dt <- hits_dt[hits_dt$isAnchor==TRUE, ]
    genome_1 <- anchors_dt$genome1[1]
    genome_2 <- anchors_dt$genome2[1]

    colnames(anchors_dt)[colnames(anchors_dt)==c("chr1")] <- paste0("chr_",genome_1)
    colnames(anchors_dt)[colnames(anchors_dt)==c("start1")] <- paste0("start_",genome_1)
    colnames(anchors_dt)[colnames(anchors_dt)==c("end1")] <- paste0("end_",genome_1)
    colnames(anchors_dt)[colnames(anchors_dt)==c("id1")] <- paste0("id_",genome_1)

    colnames(anchors_dt)[colnames(anchors_dt)==c("chr2")] <- paste0("chr_",genome_2)
    colnames(anchors_dt)[colnames(anchors_dt)==c("start2")] <- paste0("start_",genome_2)
    colnames(anchors_dt)[colnames(anchors_dt)==c("end2")] <- paste0("end_",genome_2)
    colnames(anchors_dt)[colnames(anchors_dt)==c("id2")] <- paste0("id_",genome_2)

    merge_dt <- merge(merge_dt, anchors_dt, by = intersect(names(merge_dt), names(anchors_dt)))

    which_relevant <- which(colnames(merge_dt) %in% c(paste0(c("chr_", "start_", "end_", "id_"),genome_1),
                                     paste0(c("chr_", "start_", "end_", "id_"),genome_2)))

    merge_dt <- merge_dt[, ..which_relevant]

  }

  return(merge_dt)
}



#' map_bins_to_anchors
#'
#' @param heal_list
#' @param genespace_dir
#' @param n_cores
#' @param bin_size
#'
#' @return
#'
map_bins_to_anchors <- function(heal_list, genespace_dir){

  map_dt_list <- parse_genespace_input(genespace_dir)
  syn_hits_list <- map_dt_list$syn_hits

  bins_to_map_list <- get_bins_to_map(heal_list, map_dt_list$block_coord)

  anchors_dt <- get_conserved_anchors(genespace_dir)

  progenitors <- names(heal_list)

  syn_anchors_to_bin_dt_list <- list()
  for(prog in progenitors){

    chr <- paste0("chr_",prog)
    start <- paste0("start_",prog)
    end <- paste0("end_",prog)
    id <- paste0("id_",prog)

    gnm_anchors_ranges <- GenomicRanges::GRanges(seqnames = anchors_dt[[chr]],
                                                  ranges = IRanges::IRanges(start = anchors_dt[[start]],
                                                                            end = anchors_dt[[end]]),
                                                  gene_id = anchors_dt[[id]])

    gnm_bin_ranges <- GenomicRanges::GRanges(seqnames = bins_to_map_list[[prog]]$bin_chr,
                                                 ranges = IRanges::IRanges(start = bins_to_map_list[[prog]]$bin_start,
                                                                           end = bins_to_map_list[[prog]]$bin_end))

    gnm_overlap_bins_anchors <- GenomicRanges::findOverlaps(query = gnm_bin_ranges, subject = gnm_anchors_ranges)

    dt_overlap_bins_anchors <- data.table::as.data.table(gnm_overlap_bins_anchors) #Because queryHits does not work..?

    gnm_overlap_widths <- IRanges::width(IRanges::pintersect(gnm_bin_ranges[dt_overlap_bins_anchors$queryHits], gnm_anchors_ranges[dt_overlap_bins_anchors$subjectHits]))

    gnm_overlap_dt <- data.table::data.table(
      bin_index = dt_overlap_bins_anchors$queryHits,
      anchor_index = dt_overlap_bins_anchors$subjectHits,
      bin_start = IRanges::start(gnm_bin_ranges[dt_overlap_bins_anchors$queryHits]),
      anchor_start = IRanges::start(gnm_anchors_ranges[dt_overlap_bins_anchors$subjectHits]),
      anchor_end = IRanges::end(gnm_anchors_ranges[dt_overlap_bins_anchors$subjectHits]),
      gene_id = gnm_anchors_ranges$gene_id[dt_overlap_bins_anchors$subjectHits],
      chr = bins_to_map_list[[prog]]$bin_chr[dt_overlap_bins_anchors$queryHits],
      overlap_width = gnm_overlap_widths
    )

    syn_anchors_to_bin_dt_list[[prog]] <- gnm_overlap_dt

  }

  return(syn_anchors_to_bin_dt_list)

}


#' get_cn_per_anchor_per_sample
#'
#' @param heal_list
#' @param genespace_dir
#' @param n_cores
#' @param bin_size
#'
#' @return
#'
get_cn_per_anchor_per_sample <- function(heal_list, genespace_dir, n_cores){

  anchors_and_bins <- map_bins_to_anchors(heal_list, genespace_dir)

  sample_names <- unlist(lapply(heal_list, function(prog){setdiff(colnames(prog$bins),c("chr", "start", "end", "mappability", "gc_content"))}))
  polyploid_samples <- names(table(sample_names))[table(sample_names)==2]

  progenitors <- names(anchors_and_bins)

  anchors_dt <- get_conserved_anchors(genespace_dir)

  anchors_bins_and_cn_list <- foreach::foreach(prog=progenitors)%do%{

    anchors_dt <- get_conserved_anchors(genespace_dir)
    doParallel::registerDoParallel(n_cores)
    cn_per_sample_list <- foreach::foreach(smp=polyploid_samples)%dopar%{

      return(heal_list[[prog]]$CN[[smp]][anchors_and_bins[[prog]]$bin_index])

    }
    doParallel::stopImplicitCluster()
    cn_dt <- as.data.table(cn_per_sample_list)
    colnames(cn_dt) <- polyploid_samples

    cn_anchors_and_bins_dt <- cbind(anchors_and_bins[[prog]],cn_dt)
    return(cn_anchors_and_bins_dt)
  }

  names(anchors_bins_and_cn_list) <- progenitors
  return(anchors_bins_and_cn_list)

}


#' get_cn_alignment_by_anchors
#'
#' @param heal_list
#' @param genespace_dir
#' @param n_cores
#'
#' @return
#'
get_cn_alignment_by_anchors <- function(heal_list, genespace_dir, n_cores){

  anchors_dt <- get_conserved_anchors(genespace_dir)
  cn_per_anchor_per_sample_dt <- get_cn_per_anchor_per_sample(heal_list, genespace_dir, n_cores)

  sample_names <- unlist(lapply(heal_list, function(prog){setdiff(colnames(prog$bins),c("chr", "start", "end", "mappability", "gc_content"))}))
  polyploid_samples <- names(table(sample_names))[table(sample_names)==2]
  progenitors <- names(heal_list)

  cn_alignment_list <- foreach::foreach(smp=polyploid_samples)%do%{

    merge_dt <- anchors_dt

    foreach::foreach(prog=progenitors)%do%{

      id <- cn_per_anchor_per_sample_dt[[prog]]$gene_id
      cn <- cn_per_anchor_per_sample_dt[[prog]][[smp]]
      cn_dt <- data.table(id, cn)
      colnames(cn_dt) <- c(paste0("id_", prog), paste0("cn_", prog))

      merge_dt <- merge(merge_dt, cn_dt, by = intersect(names(merge_dt), names(cn_dt)))

    }

    return(merge_dt)

  }

  names(cn_alignment_list) <- polyploid_samples
  return(cn_alignment_list)

}




#' Title
#'
#' @param heal_list
#' @param tmp_map_dt
#' @param ref_gnm
#' @param ref_chr
#' @param alt_gnm
#' @param alt_chr
#' @param bin_size
#'
#' @return
#' @export
#'
#' @importFrom data.table :=
dtw_na <- function(heal_list, tmp_map_dt, ref_gnm, ref_chr, alt_gnm, alt_chr, bin_size){

  smp_medians <- get_sample_stats(heal_list)

  alt_existing_bins <- heal_list[[alt_gnm]]$CN$start[heal_list[[alt_gnm]]$CN$chr==alt_chr]

  sample_names <- unlist(lapply(heal_list, function(prog){setdiff(colnames(prog$bins),c("chr", "start", "end", "mappability", "gc_content"))}))
  polyploid_samples <- names(table(sample_names))[table(sample_names)==2]

  for(smp in polyploid_samples){
    tmp_map_dt[,paste0("aligned_",smp) := rep("anchors",nrow(tmp_map_dt))]
  }

  indices <- 1:nrow(tmp_map_dt)
  na_indices <- indices[tmp_map_dt$alt_bin=="REF_ANCHOR_NOT_IN_BIN" | tmp_map_dt$alt_bin=="ALT_ANCHOR_NOT_IN_BIN"]
  if (length(na_indices) == 0) {
    map_dt <- data.table::data.table(ref_bin=tmp_map_dt$ref_bin)
    map_dt$ref_chr <- rep(ref_chr, nrow(map_dt))
    map_dt$alt_chr <- rep(alt_chr, nrow(map_dt))
    for(smp in polyploid_samples){
      map_dt[[smp]] <- tmp_map_dt$alt_bin
      smp_col_name <- paste0("aligned_",smp)
      col_vector <- rep("anchors",nrow(tmp_map_dt))
      data.table::setDT(map_dt)
      data.table::set(map_dt, j = smp_col_name, value = col_vector)
    }

    return(map_dt)
  }

  no_na_indices <- indices[tmp_map_dt$alt_bin!="REF_ANCHOR_NOT_IN_BIN" & tmp_map_dt$alt_bin!="ALT_ANCHOR_NOT_IN_BIN"]

  # If start and end of block has missing bin, truncate block to first available
  if(min(na_indices)==1){
    first_not_na <- min(no_na_indices)

    if(max(na_indices)==nrow(tmp_map_dt)){
      last_not_na <- max(no_na_indices)
      tmp_map_dt <- tmp_map_dt[first_not_na:last_not_na,]
    }else{
      last_not_na <- nrow(tmp_map_dt)
      tmp_map_dt <- tmp_map_dt[first_not_na:last_not_na,]
    }
  }else if(max(na_indices)==nrow(tmp_map_dt)){
    first_not_na <- 1
    last_not_na <- max(no_na_indices)
    tmp_map_dt <- tmp_map_dt[first_not_na:last_not_na,]
  }

  # recompute na indices
  indices <- 1:nrow(tmp_map_dt)
  na_indices <- indices[tmp_map_dt$alt_bin=="REF_ANCHOR_NOT_IN_BIN" | tmp_map_dt$alt_bin=="ALT_ANCHOR_NOT_IN_BIN"]

  if (length(na_indices) == 0) {

    map_dt <- data.table::data.table(ref_bin=tmp_map_dt$ref_bin)
    map_dt$ref_chr <- rep(ref_chr, nrow(map_dt))
    map_dt$alt_chr <- rep(alt_chr, nrow(map_dt))
    for(smp in polyploid_samples){
      map_dt[[smp]] <- tmp_map_dt$alt_bin
      smp_col_name <- paste0("aligned_",smp)
      col_vector <- rep("anchors",nrow(tmp_map_dt))
      data.table::setDT(map_dt)
      data.table::set(map_dt, j = smp_col_name, value = col_vector)
    }

    return(map_dt)
  }

  # Identify the start and end of each NA sequence
  na_runs <- split(na_indices, cumsum(c(1, diff(na_indices) != 1)))

  per_run_replacement_lists <- lapply(na_runs, function(na_run){
    print(na_run)
    start_na <- na_run[1]
    end_na <- na_run[length(na_run)]

    ref_available_bins <- tmp_map_dt$ref_bin[start_na:end_na]

    alt_prev_position <- as.numeric(tmp_map_dt$alt_bin[start_na - 1])
    alt_next_position <- as.numeric(tmp_map_dt$alt_bin[end_na + 1])
    if((alt_next_position-alt_prev_position)>0){
      alt_available_bins <- alt_existing_bins[alt_existing_bins>=alt_prev_position & alt_existing_bins<=alt_next_position]
    }else{
      alt_available_bins <- alt_existing_bins[alt_existing_bins<=alt_prev_position & alt_existing_bins>=alt_next_position]
    }

    if(length(alt_available_bins)==1){

      replacement_list <- foreach::foreach(smp=polyploid_samples)%do%{
        smp_col_name <- paste0("aligned_",smp)
        tmp_map_dt[start_na:end_na , (smp_col_name) := rep("single_alt", length(start_na:end_na))]
        replacement_values <- rep(alt_available_bins, length(na_run))
        return(replacement_values)

      }
      names(replacement_list) <- polyploid_samples
      return(replacement_list)

    }else{

      ref_cn_at_bins <-  dplyr::filter(heal_list[[ref_gnm]]$CN, start %in% ref_available_bins & chr %in% ref_chr)
      alt_cn_at_bins <-  dplyr::filter(heal_list[[alt_gnm]]$CN, start %in% alt_available_bins & chr %in% alt_chr)
      ref_count_at_bins <-  dplyr::filter(heal_list[[ref_gnm]]$bins, start %in% ref_available_bins & chr %in% ref_chr)
      alt_count_at_bins <-  dplyr::filter(heal_list[[alt_gnm]]$bins, start %in% alt_available_bins & chr %in% alt_chr)

      replacement_list <- foreach::foreach(smp=polyploid_samples)%do%{

        smp_col_name <- paste0("aligned_",smp)
        ref_cn_vec <- ref_cn_at_bins[[smp]]
        alt_cn_vec <- alt_cn_at_bins[[smp]]

        # normalize count values
        ref_count_vec <- abs(ref_count_at_bins[[smp]]/smp_medians[[smp]]-1)
        alt_count_vec <- abs(alt_count_at_bins[[smp]]/smp_medians[[smp]]-1)

        # Replace NAs (count outliers) by mean of the non-na values (inspired by CBS method)
        #cat("Why not any number? and why x2?")
        #cat("why do I always get counts? Sometimes not needed?")
        if(sum(is.na(ref_count_vec))==length(ref_count_vec)){
          if(sum(is.na(alt_count_vec)==length(alt_count_vec))){

            ref_count_vec <- rep(smp_medians[[smp]]*2, length(ref_count_vec))
            alt_count_vec <- rep(smp_medians[[smp]]*2, length(alt_count_vec))
          }else{
            ref_count_vec <- rep(mean(na.omit(alt_count_vec), length(ref_count_vec)))
          }

        }else if(sum(is.na(alt_count_vec))==length(alt_count_vec)){
          alt_count_vec <- rep(mean(na.omit(ref_count_vec), length(alt_count_vec)))
        }
        ref_count_vec[is.na(ref_count_vec)] <- mean(na.omit(ref_count_vec))
        alt_count_vec[is.na(alt_count_vec)] <- mean(na.omit(alt_count_vec))

        # If one of the CN is uninformative i.e. constant
        if(length(unique(ref_cn_vec))==1 | length(unique(alt_cn_vec))==1){
          # If there is discordance i.e. only 1 is constant or both constant are not according
          if(sum(abs(unique(ref_cn_vec)-2)!=abs(unique(alt_cn_vec)-2))>0){

            tmp_map_dt[start_na:end_na ,  (smp_col_name) := rep("dtw_counts", length(start_na:end_na))]

            # Align based on normalize counts

            count_alignment <- dtw::dtw(ref_count_vec, alt_count_vec)

            # For each ref bin mapping multiple alt chose single best match
            ref_section <- split(count_alignment$index1, cumsum(c(1, diff(count_alignment$index1) != 0)))
            alt_section <- split(count_alignment$index2, cumsum(c(1, diff(count_alignment$index1) != 0)))


            replacement_values <- c()
            for(i in 1:length(ref_section)){
              if(length(ref_count_vec[ref_section[[i]]])>1){
                which_less <- which.min(abs(ref_count_vec[ref_section[[i]]]-alt_count_vec[alt_section[[i]]]))
                best_alt <- alt_section[[i]][which_less]
                replacement_values <- c(replacement_values, alt_available_bins[best_alt])

              }else{
                replacement_values <- c(replacement_values, alt_available_bins[alt_section[[i]]])
              }
            }
            return(replacement_values)

          # Concordance i.e. unique and concordant
          }else{

            tmp_map_dt[start_na:end_na , (smp_col_name) := rep("spread_concordant_CN", length(start_na:end_na))]

            if(length(ref_available_bins)==(length(alt_available_bins)-2)){
              replacement_values <- alt_available_bins[2:(length(alt_available_bins)-1)]

            }else if(length(ref_available_bins)>(length(alt_available_bins)-2)){

              if(length(ref_available_bins)==(length(alt_available_bins)-1)){
                replacement_values <- alt_available_bins[1:length(ref_available_bins)]
                return(replacement_values)

              }else{
                long_vec_len <- length(ref_available_bins)
                short_vec_len <- length(alt_available_bins)
                indices <- round(seq(1, long_vec_len, length.out = short_vec_len))
                replacement_values <- rep(NA, long_vec_len)
                replacement_values[indices] <- alt_available_bins
                replacement_values <- zoo::na.locf(replacement_values)
                return(replacement_values)
              }

            }else if(length(ref_available_bins)<(length(alt_available_bins)-2)){

              indices <- seq(1, (length(alt_available_bins)-2), length.out=length(ref_available_bins))
              replacement_values <- alt_available_bins[2:(length(alt_available_bins)-1)][indices]
              return(replacement_values)
            }
          }

        # Informative CN
        }else{

          tmp_map_dt[start_na:end_na , (smp_col_name) := rep("dtw_cn", length(start_na:end_na))]

          cn_alignment <- dtw::dtw(ref_cn_vec, alt_cn_vec)

          ref_section <- split(cn_alignment$index1, cumsum(c(1, diff(cn_alignment$index1) != 0)))
          alt_section <- split(cn_alignment$index2, cumsum(c(1, diff(cn_alignment$index1) != 0)))

          replacement_values <- c()
          for(i in 1:length(ref_section)){
            if(length(ref_count_vec[ref_section[[i]]])>1){

              which_less <- which.min(abs(ref_count_vec[ref_section[[i]]]-alt_count_vec[alt_section[[i]]]))
              best_alt <- alt_section[[i]][which_less]
              replacement_values <- c(replacement_values, alt_available_bins[best_alt])

            }else{
              replacement_values <- c(replacement_values, alt_available_bins[alt_section[[i]]])
            }
          }
          return(replacement_values)

        }

      }

      names(replacement_list) <- polyploid_samples
      return(replacement_list)
    }
  })

  replacement_df <- do.call(rbind, lapply(per_run_replacement_lists,as.data.frame))

  map_dt <- data.table::data.table(ref_bin=tmp_map_dt$ref_bin)
  map_dt$ref_chr <- rep(ref_chr, nrow(map_dt))
  map_dt$alt_chr <- rep(alt_chr, nrow(map_dt))

  for(smp in polyploid_samples){
    smp_col_name <- paste0("aligned_",smp)
    col_vector <- tmp_map_dt[[smp_col_name]]
    data.table::setDT(map_dt)
    data.table::set(map_dt, j = smp, value = suppressWarnings(as.numeric(tmp_map_dt$alt_bin)))
    #map_dt[[smp]] <- suppressWarnings(as.numeric(tmp_map_dt$alt_bin))
    data.table::set(map_dt, j = smp_col_name, value = col_vector)
  }

  if(length(unlist(na_runs))!=nrow(replacement_df)){
    cat("ERROR: Replacement and NA not same length.")
    return()
  }

  for(i in 1:nrow(replacement_df)){
    position <- unlist(na_runs)[i]

    for(smp in polyploid_samples){
      bin_to_map <- replacement_df[i,smp]
      map_dt[position, (smp) := bin_to_map]
      }
  }

  return(map_dt)
}


align_blocks <- function(blk_dt, heal_list, ref_gnm, alt_gnm, ref_anchors_dt, alt_anchors_dt, bin_size, n_cores){
  cat("If constant and discordant, do we need dtw?")
  doParallel::registerDoParallel(n_cores)

  map_per_blk_list <- foreach::foreach(i=(1:nrow(blk_dt)))%dopar%{

    blk_id <- blk_dt$blkID[i]
    cat(paste0("Processing syntenic block ",blk_id,".\n"))
    ref_chr <- blk_dt$chr1[i]
    ref_start <- blk_dt$startBp1[i]
    ref_end <- blk_dt$endBp1[i]

    alt_chr <- blk_dt$chr2[i]
    alt_start <- blk_dt$minBp2[i]
    alt_end <- blk_dt$maxBp2[i]

    #############################
    # Get bins within the block.#
    #############################

    ### Ref ###
    dt_ref <- heal_list[[ref_gnm]]$bins[heal_list[[ref_gnm]]$bins$chr==ref_chr,]
    ref_bin_centers <- dt_ref$start+(bin_size/2)
    ref_blk_orientation <- sign(ref_end-ref_start)

    if(ref_blk_orientation==1){
      n_bins_in_blk <- sum(ref_bin_centers>ref_start & ref_bin_centers<ref_end)
      if(n_bins_in_blk==0){
        return()
      }
    }else if(ref_blk_orientation==-1){
      n_bins_in_blk <- sum(ref_bin_centers<ref_start & ref_bin_centers>ref_end)
      if(n_bins_in_blk==0){
        return()
      }
    }else{
      cat(paste("ERROR: blk borders for",blk_id,"not valid. Exiting.."))
      return()
    }

    # Gene center within a bin
    if(min(abs(ref_bin_centers-ref_start))<(bin_size/2)){
      ref_blk_bin_start <- dt_ref$start[which.min(abs(ref_bin_centers-ref_start))]

      # If not, start is the first available bin towards center of the block.
    }else{
      distance_to_strt_vec <- ref_bin_centers-ref_start
      if(ref_blk_orientation==1){
        only_positive_dist <- distance_to_strt_vec[distance_to_strt_vec>1]
        min_positive <- which.min(only_positive_dist) # this should always be one.
        ref_blk_bin_start <- dt_ref$start[distance_to_strt_vec==only_positive_dist[min_positive]]

      }else if(ref_blk_orientation==-1){
        only_negative_dist <- distance_to_strt_vec[distance_to_strt_vec<1]
        max_negative <- which.max(only_negative_dist) # this should always be length only negative.
        ref_blk_bin_start <- dt_ref$start[distance_to_strt_vec==only_negative_dist[max_negative]]

      }else{
        cat(paste("ERROR: blk borders for",blk_id,"not valid. Exiting.."))
        return()
      }
    }


    # Gene center within a bin
    if(min(abs(ref_bin_centers-ref_end))<(bin_size/2)){
      ref_blk_bin_end <- dt_ref$start[which.min(abs(ref_bin_centers-ref_end))]

      # If not, end is the first available bin towards center of the block.
    }else{
      distance_to_end_vec <- ref_bin_centers-ref_end
      if(ref_blk_orientation==-1){
        only_positive_dist <- distance_to_end_vec[distance_to_end_vec>1]
        min_positive <- which.min(only_positive_dist) # this should always be one.
        ref_blk_bin_end <- dt_ref$start[distance_to_end_vec==only_positive_dist[min_positive]]

      }else if(ref_blk_orientation==1){
        only_negative_dist <- distance_to_end_vec[distance_to_end_vec<1]
        max_negative <- which.max(only_negative_dist) # this should always be length only negative.
        ref_blk_bin_end <- dt_ref$start[distance_to_end_vec==only_negative_dist[max_negative]]

      }
    }

    if(ref_blk_orientation==1){
      ref_start_vec <- dt_ref$start[dt_ref$start>=ref_blk_bin_start & dt_ref$start<=ref_blk_bin_end]
    }else{
      ref_start_vec <- dt_ref$start[dt_ref$start<=ref_blk_bin_start & dt_ref$start>=ref_blk_bin_end]
    }
    ref_end_vec <- ref_start_vec+bin_size
    ref_center_vec <- ref_start_vec+(bin_size/2)
    ref_chromo_vec <- rep(ref_chr,length(ref_start_vec))


    ### Alt ###
    dt_alt <- heal_list[[alt_gnm]]$bins[heal_list[[alt_gnm]]$bins$chr==alt_chr,]
    alt_bin_centers <- dt_alt$start+(bin_size/2)
    alt_blk_orientation <- sign(alt_end-alt_start)

    if(alt_blk_orientation==1){
      n_bins_in_blk <- sum(alt_bin_centers>alt_start & alt_bin_centers<alt_end)
      if(n_bins_in_blk==0){
        return()
      }
    }else if(alt_blk_orientation==-1){
      n_bins_in_blk <- sum(alt_bin_centers<alt_start & alt_bin_centers>alt_end)
      if(n_bins_in_blk==0){
        return()
      }
    }else{
      cat(paste("ERROR: blk borders for",blk_id,"not valid. Exiting.."))
      return()
    }

    # Gene center within a bin
    if(min(abs(alt_bin_centers-alt_start))<(bin_size/2)){
      alt_blk_bin_start <- dt_alt$start[which.min(abs(alt_bin_centers-alt_start))]

      # If not, start is the first available bin towards center of the block.
    }else{
      distance_to_strt_vec <- alt_bin_centers-alt_start
      if(alt_blk_orientation==1){
        only_positive_dist <- distance_to_strt_vec[distance_to_strt_vec>1]
        min_positive <- which.min(only_positive_dist) # this should always be one.

        alt_blk_bin_start <- dt_alt$start[distance_to_strt_vec==only_positive_dist[min_positive]]

      }else if(alt_blk_orientation==-1){
        only_negative_dist <- distance_to_strt_vec[distance_to_strt_vec<1]
        max_negative <- which.max(only_negative_dist) # this should always be length only negative.

        alt_blk_bin_start <- dt_alt$start[distance_to_strt_vec==only_negative_dist[max_negative]]

      }else{
        cat(paste("ERROR: blk borders for",blk_id,"not valid. Exiting.."))
        return()
      }
    }


    # Gene center within a bin
    if(min(abs(alt_bin_centers-alt_end))<(bin_size/2)){
      alt_blk_bin_end <- dt_alt$start[which.min(abs(alt_bin_centers-alt_end))]

      # If not, end is the first available bin towards center of the block.
    }else{
      distance_to_end_vec <- alt_bin_centers-alt_end
      if(alt_blk_orientation==-1){
        only_positive_dist <- distance_to_end_vec[distance_to_end_vec>1]
        min_positive <- which.min(only_positive_dist) # this should always be one.

        alt_blk_bin_end <- dt_alt$start[distance_to_end_vec==only_positive_dist[min_positive]]

      }else if(alt_blk_orientation==1){
        only_negative_dist <- distance_to_end_vec[distance_to_end_vec<1]
        max_negative <- which.max(only_negative_dist) # this should always be length only negative.

        alt_blk_bin_end <- dt_alt$start[distance_to_end_vec==only_negative_dist[max_negative]]

      }
    }


    if(alt_blk_orientation==1){
      alt_start_vec <- dt_alt$start[dt_alt$start>=alt_blk_bin_start & dt_alt$start<=alt_blk_bin_end]
    }else{
      alt_start_vec <- dt_alt$start[dt_alt$start<=alt_blk_bin_start & dt_alt$start>=alt_blk_bin_end]
    }

    alt_end_vec <- alt_start_vec+bin_size
    alt_center_vec <- alt_start_vec+(bin_size/2)
    alt_chromo_vec <- rep(alt_chr,length(alt_start_vec))


    ###########################
    # Align within block bins #
    ###########################
    # we use anchors for that
    cat("Maybe don't need the vectors. Just align bins based on anchors and that's it?")

    ref_anchors_in_blk <- ref_anchors_dt[ref_anchors_dt$blkID==blk_id,]
    ref_start_end <- colnames(ref_anchors_in_blk)[grepl("^(start|end)",colnames(ref_anchors_in_blk))]
    ref_anchors_in_blk$center_point <- rowMeans(ref_anchors_in_blk[,..ref_start_end])

    ref_anchors_in_blk_list <- split(ref_anchors_in_blk, seq(nrow(ref_anchors_in_blk)))

    # Assign a bin to each anchor gene (if possible)
    ref_which_bin <- unlist(lapply(ref_anchors_in_blk_list, function(row){

      # Gene center within a bin

      if(min(abs(ref_center_vec-row$center_point))==-Inf | min(abs(ref_center_vec-row$center_point))==Inf){
        while(TRUE==TRUE){
          cat(blk_id)
          print("ERROROROROR")
          }
        return()
      }else if(is.na(min(abs(ref_center_vec-row$center_point)))){
      while(TRUE==TRUE){
        cat(blk_id)
        print("ERROROROROR")
      }
        cat(blk_id)
        print("EROROROR")
        return()
      }
      if(min(abs(ref_center_vec-row$center_point))<(bin_size/2)){
        bin_start <- ref_start_vec[which.min(abs(ref_center_vec-row$center_point))]
        return(bin_start)

        # We assume no gene spans more than 2 bins such that
        # if match the following 2 ifs, these are the only bins possible.
        print('WRONG ASSUMPTION!')
        # Start within a bin
      }else if(sum(row$start1>ref_start_vec & row$start1<ref_end_vec)==1){
        bin_start <- ref_start_vec[which(row$start1>ref_start_vec & row$start1<ref_end_vec)]
        return(bin_start)

        # End within a bin
      }else if(sum(row$end1>ref_start_vec & row$end1<ref_end_vec)==1){
        bin_start <- ref_start_vec[which(row$end1>ref_start_vec & row$end1<ref_end_vec)]
        return(bin_start)

      }else{
        return("NO_BIN")
      }
    }))


    alt_anchors_in_blk <- alt_anchors_dt[alt_anchors_dt$blkID==blk_id,]
    alt_start_name <- colnames(alt_anchors_in_blk)[grepl("^start",colnames(alt_anchors_in_blk))]
    alt_anchor_start<- alt_anchors_in_blk[[alt_start_name]]
    alt_end_name <- colnames(alt_anchors_in_blk)[grepl("^end",colnames(alt_anchors_in_blk))]
    alt_anchor_end <- alt_anchors_in_blk[[alt_end_name]]
    poz_col <- c(alt_start_name,alt_end_name)
    alt_anchor_centers <- rowMeans(alt_anchors_in_blk[,..poz_col])
    names(alt_anchor_centers) <- names(alt_anchor_start) <- names(alt_anchor_end) <- ref_which_bin


    ref_to_alt_list <- sapply(ref_start_vec,function(bin){
      ref_start <- bin
      alt_mean <- mean(alt_anchor_centers[as.character(bin)])
      alt_min <- min(alt_anchor_start[as.character(bin)])
      alt_max <- max(alt_anchor_end[as.character(bin)])
      list(ref_start=ref_start, alt_mean=alt_mean, alt_min=alt_min, alt_max=alt_max)}, simplify = FALSE)

    which_bin_alt <- unlist(lapply(ref_to_alt_list ,function(row){

      if(is.na(row$alt_mean)){
        return("REF_ANCHOR_NOT_IN_BIN")
      }else{

        # Homoeolog edges containing bins?
        alt_bins_included <- row$alt_min<alt_end_vec & row$alt_max>alt_start_vec
        if(sum(alt_bins_included)>0){ # Bins in the range
          bin_start <- alt_start_vec[which.min(abs(alt_center_vec-row$alt_mean))]
          return(bin_start) # return bin closest to center

        }else{
          return("ALT_ANCHOR_NOT_IN_BIN")
        }
      }
    }))

    tmp_map_dt <- data.table::data.table(ref_bin=ref_start_vec,alt_bin=which_bin_alt)

    map_dt <- dtw_na(tmp_map_dt = tmp_map_dt, heal_list = heal_list, ref_gnm = ref_gnm,
                    alt_gnm = alt_gnm, ref_chr = ref_chr, alt_chr = alt_chr, bin_size = bin_size )

    return(map_dt)
  }
  doParallel::stopImplicitCluster()

  return(map_per_blk_list)

}

# related to ":=" as in: https://github.com/tidyverse/dtplyr/issues/426
.datatable.aware <- TRUE


make_aln_map <- function(heal_list, genespace_dir, bin_size, n_cores){

  cn_exist <- sum(names(heal_list[[1]])=="CN")!=0
  if(cn_exist!=TRUE){
    cat("ERROR: no CN data. Exiting...")
    return()
  }

  map_dt_list <- parse_genespace_input(genespace_dir)

  samples_list <- strsplit(names(map_dt_list$syn_anchors),"_vs_")
  samples_list <- lapply(samples_list,sort)

  genomes <- unique(map_dt_list$block_coord$genome1)
  if(sum(sort(genomes)!=sort(names(heal_list)))!=0){
    cat("ERROR: Genome names in GENESPACE output and HEAL input directory do not match. Exiting..")
    return()
  }


  key_per_genome <- foreach::foreach(ref_gnm=genomes)%do%{

    other_gnms <- setdiff(genomes,ref_gnm)

    map_to_each_other <- foreach::foreach(alt_gnm=other_gnms)%do%{

      smpl_pair <- sort(c(alt_gnm, ref_gnm))

      which_anchors <- unlist(lapply(samples_list,function(vec){
        sum(vec!=smpl_pair)==0
      }))

      anchors_dt <- map_dt_list$syn_anchors[[which_anchors]]
      cat("not sure about who is anchor")
      anchors_dt <- anchors_dt[anchors_dt$isAnchor==TRUE,]

      gnm_1_cols <- colnames(anchors_dt)[grepl("1$",colnames(anchors_dt))]
      gnm_2_cols <- colnames(anchors_dt)[grepl("2$",colnames(anchors_dt))]
      gnm_1_keep_cols <- setdiff(colnames(anchors_dt), gnm_2_cols)
      gnm_2_keep_cols <- setdiff(colnames(anchors_dt), gnm_1_cols)

      if(anchors_dt$genome1[1]==ref_gnm){
        ref_anchors_dt <- anchors_dt[, ..gnm_1_keep_cols]
        alt_anchors_dt <- anchors_dt[, ..gnm_2_keep_cols]

      }else if(anchors_dt$genome1[1]==alt_gnm){
        ref_anchors_dt <- anchors_dt[, ..gnm_2_keep_cols]
        alt_anchors_dt <- anchors_dt[, ..gnm_1_keep_cols]

      }else{
        cat("ERROR: syntenic anchors genome names don't match sample names.")
        return()
      }

      blk_dt <- map_dt_list$block_coord[map_dt_list$block_coord$genome1==ref_gnm & map_dt_list$block_coord$genome2==alt_gnm,]

      map_per_blk_list <- align_blocks(blk_dt = blk_dt, heal_list = heal_list, ref_gnm = ref_gnm, alt_gnm = alt_gnm, ref_anchors_dt = ref_anchors_dt, alt_anchors_dt = alt_anchors_dt, n_cores = n_cores, bin_size = bin_size)

      map <- data.table::rbindlist(map_per_blk_list)

      map$ref_bin <- as.numeric(map$ref_bin)

      sample_names <- unlist(lapply(heal_list, function(prog){setdiff(colnames(prog$bins),c("chr", "start", "end", "mappability", "gc_content"))}))
      polyploid_samples <- names(table(sample_names))[table(sample_names)==2]
      for(smp in polyploid_samples){
        map[[smp]] <- as.numeric(map[[smp]])
      }

      return(map)
    }
  names(map_to_each_other) <- other_gnms
  return(map_to_each_other)
  }
  names(key_per_genome) <- genomes
  return(key_per_genome)
}


per_sample_alignment <- function(heal_list, aln_map, n_cores=1, prog_ploidy=2){

  cn_exist <- sum(names(heal_list[[1]])=="CN")!=0
  if(cn_exist!=TRUE){
    cat("ERROR: no CN data. Exiting...")
    return()
  }

  sample_names <- unlist(lapply(heal_list, function(prog){setdiff(colnames(prog$bins),c("chr", "start", "end", "mappability", "gc_content"))}))
  polyploid_samples <- names(table(sample_names))[table(sample_names)==2]

  subgenomes <- names(aln_map)

  doParallel::registerDoParallel(n_cores)
  per_sample_list <- foreach::foreach(smp=polyploid_samples)%dopar%{

    ref_dt_list <- foreach::foreach(ref=subgenomes)%do%{

      alt_gnms <- setdiff(subgenomes, ref)

      col_ref_chr_name <- paste0(ref, "_chr")
      col_ref_pos_name <- paste0(ref, "_start")
      col_ref_count_name <- paste0(ref, "_count")
      col_ref_cn_name <- paste0(ref, "_cn")
      col_ref_gc_name <- paste0(ref, "_gc")
      col_ref_mapa_name <- paste0(ref, "_mappability")

      alt_dts <- foreach::foreach(alt=alt_gnms)%do%{

        map <- aln_map[[ref]][[alt]]

        ref_cn_vec <- c()
        ref_count_vec <- c()
        ref_gc_vec <- c()
        ref_mapa_vec <- c()

        alt_cn_vec <- c()
        alt_count_vec <- c()
        alt_gc_vec <- c()
        alt_mapa_vec <- c()

        for(i in 1:nrow(map)){

          ref_chr <- map$ref_chr[i]
          alt_chr <- map$alt_chr[i]
          ref_pos <- map$ref_bin[i]
          alt_pos <- map[[smp]][i]

          ref_cn <- heal_list[[ref]]$CN[[smp]][heal_list[[ref]]$CN$chr==ref_chr & heal_list[[ref]]$CN$start==ref_pos]
          alt_cn <- heal_list[[ref]]$CN[[smp]][heal_list[[ref]]$CN$chr==ref_chr & heal_list[[ref]]$CN$start==ref_pos]

          ref_count <- heal_list[[alt]]$bins[[smp]][heal_list[[alt]]$bins$chr==alt_chr & heal_list[[alt]]$bins$start==alt_pos]
          alt_count <- heal_list[[alt]]$bins[[smp]][heal_list[[alt]]$bins$chr==alt_chr & heal_list[[alt]]$bins$start==alt_pos]

          ref_gc <- heal_list[[ref]]$bins$gc_content[heal_list[[ref]]$bins$chr==ref_chr & heal_list[[ref]]$bins$start==ref_pos]
          alt_gc <- heal_list[[ref]]$bins$gc_content[heal_list[[ref]]$bins$chr==ref_chr & heal_list[[ref]]$bins$start==ref_pos]

          ref_mapa <- heal_list[[alt]]$bins$mappability[heal_list[[alt]]$bins$chr==alt_chr & heal_list[[alt]]$bins$start==alt_pos]
          alt_mapa <- heal_list[[alt]]$bins$mappability[heal_list[[alt]]$bins$chr==alt_chr & heal_list[[alt]]$bins$start==alt_pos]

          ref_cn_vec <- c(ref_cn_vec, ref_cn)
          ref_count_vec <- c(ref_count_vec, ref_count)
          ref_gc_vec <- c(ref_gc_vec, ref_gc)
          ref_mapa_vec <- c(ref_mapa_vec, ref_mapa)

          alt_cn_vec <- c(alt_cn_vec, alt_cn)
          alt_count_vec <- c(alt_count_vec, alt_count)
          alt_gc_vec <- c(alt_gc_vec, alt_gc)
          alt_mapa_vec <- c(alt_mapa_vec, alt_mapa)

        }

        col_alt_chr_name <- paste0(alt, "_chr")
        col_alt_pos_name <- paste0(alt, "_start")
        col_alt_count_name <- paste0(alt, "_count")
        col_alt_cn_name <- paste0(alt, "_cn")
        col_alt_gc_name <- paste0(alt, "_gc")
        col_alt_mapa_name <- paste0(alt, "_mappability")

        colnames_dt <- c(col_ref_chr_name,
                         col_ref_pos_name,
                         col_ref_count_name,
                         col_ref_cn_name,
                         col_ref_gc_name,
                         col_ref_mapa_name,
                         col_alt_chr_name,
                         col_alt_pos_name,
                         col_alt_count_name,
                         col_alt_cn_name,
                         col_alt_gc_name,
                         col_alt_mapa_name)

        key_cols <- colnames_dt <- c(col_ref_chr_name,
                                     col_ref_pos_name,
                                     col_alt_chr_name,
                                     col_alt_pos_name)

        aln_dt <- data.table::as.data.table(
          setNames(list(map$ref_chr, map$ref_bin, ref_count_vec,
                        ref_cn_vec, ref_gc_vec, ref_mapa_vec,
                        map$alt_chr, map[[smp]], alt_count_vec,
                        alt_cn_vec, alt_gc_vec, alt_mapa_vec),
                   colnames_dt))
        data.table::setkeyv(aln_dt, key_cols)

        return(aln_dt)
      }

      columns_to_merge_by <- c(col_ref_chr_name, col_ref_pos_name)
      map <- Reduce(function(x, y) merge(x, y, by=columns_to_merge_by), alt_dts)

      return(map)
    }

    names(ref_dt_list) <- subgenomes
    return(ref_dt_list)
  }
  doParallel::stopImplicitCluster()

  names(per_sample_list) <- polyploid_samples
  return(per_sample_list)
}





# #
# align_genes <- function(){}
# map_it <- parse_genespace_input(genespace_dir = "/srv/kenlab/kenji/heal_dirs/test_dir/genespace/")
# map_dt_list
#
# ratios <- foreach::foreach(i=1:nrow(map_it$syn_anchors$A.halleri_vs_A.lyrata))%do%{
#   id1 <- map_it$syn_anchors$A.halleri_vs_A.lyrata$id1[i]
#   count_1 <- cn_list$A.halleri$genes$HM_RS2K_G1_2[cn_list$A.halleri$genes$GeneID==id1]
#
#   id2 <- map_it$syn_anchors$A.halleri_vs_A.lyrata$id2[i]
#   count_2 <- cn_list$A.lyrata$genes$HM_RS2K_G1_2[cn_list$A.lyrata$genes$GeneID==id2]
#   if(count_1>count_2){
#     return(count_2/count_1)
#   }else{
#     return(-(count_1/count_2))
#   }
#
# }
# ratio_vec <- unlist(ratios)
#
# hist(ratio_vec, xlim=c(-1,1),breaks = 80,main="histogram of gene ratio", col="chartreuse4")#,# ylim=c(0,3000))
# cn_list$A.lyrata$genes
# cn_list$A.halleri$genes
#


