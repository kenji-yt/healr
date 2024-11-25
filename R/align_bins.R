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

  # Pick on genome as reference and keep only data tables containing this genome.
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
map_bins_to_anchors <- function(heal_list, genespace_dir, n_cores){

  map_dt_list <- parse_genespace_input(genespace_dir)
  syn_hits_list <- map_dt_list$syn_hits

  sample_names <- unlist(lapply(heal_list, function(prog){setdiff(colnames(prog$bins),c("chr", "start", "end", "mappability", "gc_content"))}))
  polyploid_samples <- names(table(sample_names))[table(sample_names)==2]

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

    gnm_bin_ranges <- GenomicRanges::GRanges(seqnames = heal_list[[prog]]$CN$chr,
                                                 ranges = IRanges::IRanges(start = heal_list[[prog]]$CN$start,
                                                                           end = heal_list[[prog]]$CN$end))

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
      chr = as.character(GenomicRanges::seqnames(gnm_bin_ranges[dt_overlap_bins_anchors$queryHits])),
      gc_content = heal_list[[prog]]$CN$gc_content[dt_overlap_bins_anchors$queryHits],
      overlap_width = gnm_overlap_widths
    )

    doParallel::registerDoParallel(n_cores)
    cn_per_sample_list <- foreach::foreach(smp=polyploid_samples)%dopar%{

      return(heal_list[[prog]]$CN[[smp]][dt_overlap_bins_anchors$queryHits])

    }
    doParallel::stopImplicitCluster()

    cn_dt <- data.table::as.data.table(cn_per_sample_list)
    colnames(cn_dt) <- polyploid_samples

    cn_anchors_and_bins_dt <- cbind(gnm_overlap_dt,cn_dt)
    data.table::setkey(cn_anchors_and_bins_dt, chr,  bin_start, anchor_start, anchor_end, gene_id, gc_content)

    syn_anchors_to_bin_dt_list[[prog]] <- cn_anchors_and_bins_dt

  }

  return(syn_anchors_to_bin_dt_list)

}


#' Create data table with anchors and copy number at the anchors
#'
#' @param heal_list
#' @param genespace_dir
#' @param n_cores
#'
#' @return
#'
get_cn_alignment_by_anchors <- function(heal_list, genespace_dir, n_cores, prog_ploidy=2){

  cn_exist <- sum(names(heal_list[[1]])=="CN")!=0
  if(cn_exist!=TRUE){
    cat("ERROR: no CN data. Exiting...")
    return()
  }

  anchors_dt <- get_conserved_anchors(genespace_dir)
  cn_anchors_and_bins <- map_bins_to_anchors(heal_list, genespace_dir, n_cores)

  sample_names <- unlist(lapply(heal_list, function(prog){setdiff(colnames(prog$bins),c("chr", "start", "end", "mappability", "gc_content"))}))
  polyploid_samples <- names(table(sample_names))[table(sample_names)==2]
  progenitors <- names(heal_list)

  total_ploidy <- length(progenitors)*prog_ploidy

  cat("Likely very ineficient to subset whole dt for each anchor...")
  cn_alignment_list <- foreach::foreach(smp=polyploid_samples)%dopar%{

    doParallel::registerDoParallel(n_cores)

    cn_per_anchor_pair_list <- foreach::foreach(i=1:nrow(anchors_dt))%dopar%{

      cn_at_anchor_dt_list <- foreach::foreach(prog=progenitors)%do%{

        anchor_gene <- anchors_dt[[paste0("id_",prog)]][i]

        which_rows <- cn_anchors_and_bins[[prog]]$gene_id==anchor_gene

        cn_at_anchor_dt <- cn_anchors_and_bins[[prog]][which_rows,]

        return(cn_at_anchor_dt)

      }
      names(cn_at_anchor_dt_list) <- progenitors

      # If the anchor overlaps no bin (eg. bin filtered out)
      if(sum(lapply(cn_at_anchor_dt_list, nrow)==0)!=0){
        return()

      }else{

        unique_cn <- lapply(cn_at_anchor_dt_list, function(dt){ unique(dt[[smp]])})
        bin_indexes <-lapply(cn_at_anchor_dt_list, function(dt){
          paste(dt$bin_index, collapse = ",")})
        bin_index_dt <- data.table::as.data.table(bin_indexes)
        colnames(bin_index_dt) <- paste0("bin_index_",colnames(bin_index_dt))

        # Unique CN in all genomes
        if(sum(unlist(lapply(unique_cn,length)))==length(unique_cn)){

          unique_cn_dt <- data.table::as.data.table(unique_cn)
          colnames(unique_cn_dt) <- paste0("cn_",colnames(unique_cn_dt))

          output_dt <- cbind(anchors_dt[i,], unique_cn_dt, bin_index_dt, data.table::data.table(method="unique"))
          return(output_dt)

        # Find most concordant CN combination (heuristic solution)
        }else{

          cartesian_product <- data.table::as.data.table(expand.grid(unique_cn))
          most_concordant <- which(abs(rowSums(cartesian_product)-total_ploidy)==min(abs(rowSums(cartesian_product)-total_ploidy)))

          concordant_combinations <- cartesian_product[most_concordant,]

          # if multiple equally concordant combinations, pick the one with most overlap
          # Note that we consider copy number of same value at start and end of the gene as equivalent
          if(nrow(concordant_combinations)>1){

            overlaps <- apply(concordant_combinations, 1, function(row) {

              overlaps <- foreach::foreach(prog=names(row))%do%{
                which_row <- cn_at_anchor_dt_list[[prog]][[smp]]==row[[prog]]
                overlap <- sum(cn_at_anchor_dt_list[[prog]]$overlap_width[which_row])
                return(overlap)
              }
              sum(unlist(overlaps))

            })

            unique_cn_dt <- concordant_combinations[which.max(overlaps)]
            colnames(unique_cn_dt) <- paste0("cn_",colnames(unique_cn_dt))

            output_dt <- cbind(anchors_dt[i,], unique_cn_dt, bin_index_dt, data.table::data.table(method="multiple_overlap"))
            return(output_dt)

          }else{

            unique_cn_dt <- concordant_combinations
            colnames(unique_cn_dt) <- paste0("cn_",colnames(unique_cn_dt))

            output_dt <- cbind(anchors_dt[i,], unique_cn_dt, bin_index_dt, data.table::data.table(method="multiple_unique"))
            return(output_dt)
          }

        }
      }
    }
    doParallel::stopImplicitCluster()

    alignment <- data.table::rbindlist(cn_per_anchor_pair_list)
    return(alignment)

  }

  names(cn_alignment_list) <- polyploid_samples
  return(cn_alignment_list)

}



#' Get data for DBSCAN copy number re-evaluation
#'
#' @param alignment
#' @param heal_list
#' @param n_cores
#' @param prog_ploidy
#' @param n_points The n parameter for MASS::kde2d(). "Number of grid points in each direction. Can be scalar or a length-2 integer vector". Default is 1000.
#'
#' @return
#' @export
#'
get_concordant_density <- function(alignment, heal_list, n_cores, prog_ploidy=2, n_points=1000){

  polyploid_samples <- names(alignment)
  progenitors <- names(heal_list)
  total_ploidy <- length(progenitors)*prog_ploidy

  cn_columns <- grep("cn_", colnames(alignment[[1]]))

  doParallel::registerDoParallel(n_cores)
  cn_count_type_dt_list <- foreach::foreach(smp=polyploid_samples)%dopar%{

    concordant_and_unique <- rowSums(alignment[[smp]][, ..cn_columns])==total_ploidy & alignment[[smp]]$method=="unique"

    dt_per_prog_list <- foreach::foreach(prog=progenitors)%do%{

      bin_index_col <- paste0("bin_index_", prog)
      cn_col <- paste0("cn_", prog)

      cc_and_u_bins <- apply(alignment[[smp]][concordant_and_unique, ..bin_index_col], 1, function(row){
        strsplit(row,",")
      })

      cc_and_u_rows <- as.numeric(unique(unlist(cc_and_u_bins)))

      chr_vec <- heal_list[[prog]]$CN$chr[cc_and_u_rows]
      start_vec <- heal_list[[prog]]$CN$start[cc_and_u_rows]
      cc_and_u_dt <- data.table::data.table(chr=chr_vec, start=start_vec)
      data.table::setkey(cc_and_u_dt, chr, start)
      cc_and_u_dt <- merge(heal_list[[prog]]$bins, cc_and_u_dt, by = c("chr", "start"))
      col_keep <- c("chr", "start", "gc_content", smp)
      cc_and_u_dt <- cc_and_u_dt[, ..col_keep]
      data.table::setkey(cc_and_u_dt, chr, start)
      colnames(cc_and_u_dt) <- c("chr", "start", "gc_content", paste0("count_",smp))
      cc_and_u_dt <- merge(heal_list[[prog]]$CN, cc_and_u_dt, by = c("chr", "start", "gc_content"))
      col_keep <- c("chr", "start", "gc_content", paste0("count_",smp), smp)
      cc_and_u_dt <- cc_and_u_dt[, ..col_keep]
      colnames(cc_and_u_dt) <- c("chr", "start", "gc_content", "count", "cn")

      cc_and_u_dt <- cc_and_u_dt[!is.na(cc_and_u_dt$count), ]
      return(cc_and_u_dt)
    }

    total_dt <- data.table::rbindlist(dt_per_prog_list)
    key_colname <- c("chr", "start", "gc_content", "cn")
    data.table::setkeyv(total_dt, key_colname)

    return(total_dt)
  }
  doParallel::stopImplicitCluster()

  names(cn_count_type_dt_list) <- polyploid_samples

  doParallel::registerDoParallel(n_cores)
  density_list <- foreach::foreach(smp=polyploid_samples)%do%{

    cn_groups <- unique(cn_count_type_dt_list[[smp]]$cn)

    density_by_cn_list <- foreach::foreach(cn=cn_groups)%do%{

      which_rows <- cn_count_type_dt_list[[smp]]$cn==cn

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



# discordant <- rowSums(alignment[[smp]][, ..cn_columns])!=total_ploidy
#
# bin_col_indexes <- grep("bin_index", colnames(alignment[[smp]]))
# alignment[[smp]][discordant, ..bin_col_indexes]
# dis_bins <- apply(alignment[[smp]][discordant, ..bin_col_indexes], 1, function(row){
#   return(row)
# })
#
# dis_cn <- apply(alignment[[smp]][discordant, ], 1, function(row){
#   return(row)
# })

