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
  syn_hits_dt_list <- lapply(syn_hits_list, function(synHits_txt){
    return(data.table::as.data.table(read.table(gzfile(synHits_txt[[1]]), header = TRUE, sep = "\t")))
  })

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
#' @importFrom data.table .SD
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

  merge_dt <- syn_hits_list[which_has_ref][[1]][, .SD, .SDcols = c(ref_chr, ref_start, ref_end, ref_id)]
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
    cn_dt <- data.table::as.data.table(cn_per_sample_list)
    colnames(cn_dt) <- polyploid_samples

    cn_anchors_and_bins_dt <- cbind(anchors_and_bins[[prog]],cn_dt)
    return(cn_anchors_and_bins_dt)
  }

  names(anchors_bins_and_cn_list) <- progenitors
  return(anchors_bins_and_cn_list)

}


#' Create data table with anchors and copy number at the anchors
#'
#' @param heal_list
#' @param genespace_dir
#' @param n_cores
#'
#' @return
#'
get_cn_alignment_by_anchors <- function(heal_list, genespace_dir, n_cores){

  cn_exist <- sum(names(heal_list[[1]])=="CN")!=0
  if(cn_exist!=TRUE){
    cat("ERROR: no CN data. Exiting...")
    return()
  }

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
      cn_dt <- data.table::data.table(id, cn)
      colnames(cn_dt) <- c(paste0("id_", prog), paste0("cn_", prog))

      merge_dt <- merge(merge_dt, cn_dt, by = intersect(names(merge_dt), names(cn_dt)))

    }

    return(merge_dt)

  }

  names(cn_alignment_list) <- polyploid_samples
  return(cn_alignment_list)

}
