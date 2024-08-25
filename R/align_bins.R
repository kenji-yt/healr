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


dtw_na <- function(heal_list, tmp_map_dt, ref_name, ref_chr, alt_name, alt_chr, bin_size){

  smp_medians <- get_sample_stats(heal_list)

  cn_exist <- sum(names(heal_list[[alt_name]])=="CN")!=0
  if(cn_exist==TRUE){
    alt_existing_bins <- heal_list[[alt_name]]$CN$start[heal_list[[alt_name]]$CN$chr==alt_chr]
  }else{
    cat("ERROR: no CN data. Exiting...")
    quit()
  }

  get_sample_stats(heal_list)
  sample_names <- unlist(lapply(heal_list, function(prog){setdiff(colnames(prog$bins),c("chr", "start", "end", "mappability", "gc_content"))}))
  polyploid_samples <- names(table(sample_names))[table(sample_names)==2]


  na_indices <- which(tmp_map_dt$alt_bin=="REF_ANCHOR_NOT_IN_BIN" | tmp_map_dt$alt_bin=="ALT_ANCHOR_NOT_IN_BIN")

  if (length(na_indices) == 0) {
    return(tmp_map_dt)
  }

  # Identify the start and end of each NA sequence
  na_runs <- split(na_indices, cumsum(c(1, diff(na_indices) != 1)))

  per_run_replacement_lists <- lapply(na_runs, function(na_run){

    start_na <- na_run[1]
    end_na <- na_run[length(na_run)]

    ref_available_bins <- tmp_map_dt$ref_bin[start_na:end_na]

    alt_prev_position <- as.numeric(tmp_map_dt$alt_bin[start_na - 1])
    alt_next_position <- as.numeric(tmp_map_dt$alt_bin[end_na + 1])
    alt_available_bins <- alt_existing_bins[alt_existing_bins>=alt_prev_position & alt_existing_bins<=alt_next_position]

    if(length(alt_available_bins)==1){

      replacement_list <- foreach::foreach(smp=polyploid_samples)%do%{
        replacement_values <- rep(alt_available_bins, length(na_run))
        return(replacement_values)
      }
      names(replacement_list) <- polyploid_samples
      return(replacement_list)

    }else{

      ref_cn_at_bins <-  dplyr::filter(heal_list[[ref_name]]$CN, start %in% ref_available_bins & chr %in% ref_chr)
      alt_cn_at_bins <-  dplyr::filter(heal_list[[alt_name]]$CN, start %in% alt_available_bins & chr %in% alt_chr)
      ref_count_at_bins <-  dplyr::filter(heal_list[[ref_name]]$bins, start %in% ref_available_bins & chr %in% ref_chr)
      alt_count_at_bins <-  dplyr::filter(heal_list[[alt_name]]$bins, start %in% alt_available_bins & chr %in% alt_chr)

      replacement_list <- foreach::foreach(smp=polyploid_samples)%do%{

        ref_cn_vec <- ref_cn_at_bins[[smp]]
        alt_cn_vec <- alt_cn_at_bins[[smp]]

        ref_count_vec <- abs(ref_count_at_bins[[smp]]/smp_medians[[smp]]-1)
        alt_count_vec <- abs(alt_count_at_bins[[smp]]/smp_medians[[smp]]-1)

        # If one of the CN is uninformative i.e. constant
        if(length(unique(ref_cn_vec))==1 | length(unique(alt_cn_vec))==1){
          # If there is discordance i.e. only 1 is constant or both constant are not according
          if(sum(abs(unique(ref_cn_vec)-2)!=abs(unique(alt_cn_vec)-2))>0){
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

  cat("ERROR: No good behaviour for start and end missing matches. Go back to align bins to fix start-end situation.")
  replacement_df <- do.call(rbind, lapply(per_run_replacement_lists,as.data.frame))


  map_dt <- data.table::data.table(ref_bin=tmp_map_dt$ref_bin)
  map_dt$ref_chr <- rep(ref_chr, nrow(map_dt))
  map_dt$alt_chr <- rep(alt_chr, nrow(map_dt))

  for(smp in polyploid_samples){
    map_dt[[smp]] <- suppressWarnings(as.numeric(tmp_map_dt$alt_bin))
  }

  names_rows <- sort(as.numeric(gsub("\\.","",rownames(replacement_df))))
  na_runs_names <- (names(unlist(na_runs)))

  replacement_df$name_row <- names_rows

  for(name_r in names_rows){
    if(sum(name_r==na_runs_names)==1){
      position <- unlist(na_runs)[name_r]
      for(smp in polyploid_samples){

        map_dt[position, (smp) := replacement_df[replacement_df$name_row == name_r,smp]]

      }
    }
  }
#
#
#     for(i in 1:nrow(replacement_df)){
#       index <- replacement_df$df_position[i]
#       map_dt[[smp]][index] <- replacement_df[[smp]][i]
#     }
#   }
  return(map_dt)

}





align_bins <- function(heal_list, genespace_dir, bin_size){

  map_dt_list <- parse_genespace_input(genespace_dir)

  samples_list <- strsplit(names(map_dt_list$syn_hits),"_vs_")
  samples_list <- lapply(samples_list,sort)

  genomes <- unique(map_dt_list$block_coord$genome1)
  if(sum(sort(genomes)!=sort(names(heal_list)))!=0){
    cat("ERROR: Genome names in GENESPACE output and HEAL input directory do not match. Exiting..")
    quit()
  }


  key_per_genome <- foreach::foreach(ref_gnm=genomes)%do%{

    other_gnms <- setdiff(genomes,ref_gnm)

    map_to_each_other <- foreach::foreach(alt_gnm=other_gnms)%do%{

      smpl_pair <- sort(c(alt_gnm, ref_gnm))

      which_anchors <- unlist(lapply(samples_list,function(vec){
        sum(vec!=smpl_pair)==0
      }))

      anchors_dt <- map_dt_list$syn_hits[[which_anchors]]
      cat("not sure about who is anchor")
      anchors_dt <- anchors_dt[anchors_dt$isAnchor==TRUE,]

      gnm_1_cols <- colnames(anchors_dt)[grepl("1$",colnames(anchors_dt))]
      gnm_2_cols <- colnames(anchors_dt)[grepl("2$",colnames(anchors_dt))]
      gnm_1_keep_cols <- setdiff(colnames(anchors_dt), gnm_2_cols)
      gnm_2_keep_cols <- setdiff(colnames(anchors_dt), gnm_1_cols)

      if(anchors_dt$genome1[1]==ref_name){
        ref_anchors_dt <- anchors_dt[, ..gnm_1_keep_cols]
        alt_anchors_dt <- anchors_dt[, ..gnm_2_keep_cols]

      }else if(anchors_dt$genome1[1]==alt_name){
        ref_anchors_dt <- anchors_dt[, ..gnm_2_keep_cols]
        alt_anchors_dt <- anchors_dt[, ..gnm_1_keep_cols]

      }else{
        cat("ERROR: syntenic hits genome names don't match sample names.")
        quit()
      }

      blks <- map_dt_list$block_coord[map_dt_list$block_coord$genome1==ref_gnm & map_dt_list$block_coord$genome2==alt_gnm,]

      doParallel::registerDoParallel(n_cores)

      map_per_blk_list <- foreach::foreach(i=(1:nrow(blks)))%dopar%{
        print(i)
        blk_id <- blks$blkID[i]

        ref_chr <- blks$chr1[i]
        ref_start <- blks$startBp1[i]
        ref_end <- blks$endBp1[i]

        alt_chr <- blks$chr2[i]
        alt_start <- blks$minBp2[i]
        alt_end <- blks$maxBp2[i]

        dt_ref <- heal_list[[ref_gnm]]$bins[heal_list[[ref_gnm]]$bins$chr==ref_chr,]
        ref_bin_centers <- dt_ref$start+(bin_size/2)
        ref_bin_start <- dt_ref$start[which.min(abs(ref_bin_centers-ref_start))]
        ref_bin_end <- dt_ref$start[which.min(abs(ref_bin_centers-ref_end))]
        if(sign(ref_bin_end-ref_bin_start)==1){
          ref_start_vec <- dt_ref$start[dt_ref$start>=ref_bin_start & dt_ref$start<=ref_bin_end]
        }else{
          ref_start_vec <- dt_ref$start[dt_ref$start<=ref_bin_start & dt_ref$start>=ref_bin_end]
        }
        ref_end_vec <- ref_start_vec+bin_size
        ref_center_vec <- ref_start_vec+(bin_size/2)
        ref_chromo_vec <- rep(ref_chr,length(ref_start_vec))


        dt_alt <- heal_list[[alt_gnm]]$bins[heal_list[[alt_gnm]]$bins$chr==alt_chr,]
        alt_bin_centers <- dt_alt$start+(bin_size/2)
        alt_bin_start <- dt_alt$start[which.min(abs(alt_bin_centers-alt_start))]
        alt_bin_end <- dt_alt$start[which.min(abs(alt_bin_centers-alt_end))]
        if(sign(alt_bin_end-alt_bin_start)==1){
          alt_start_vec <- dt_alt$start[dt_alt$start>=alt_bin_start & dt_alt$start<=alt_bin_end]
        }else{
          alt_start_vec <- dt_alt$start[dt_alt$start<=alt_bin_start & dt_alt$start>=alt_bin_end]
        }
        alt_end_vec <- alt_start_vec+bin_size
        alt_center_vec <- alt_start_vec+(bin_size/2)
        alt_chromo_vec <- rep(alt_chr,length(alt_start_vec))

        # Align alt_x_vec to ref_x_vec
        # we use anchors for that
        ref_anchors_in_blk <- ref_anchors_dt[ref_anchors_dt$blkID==blk_id,]
        ref_start_end <- colnames(ref_anchors_in_blk)[grepl("^(start|end)",colnames(ref_anchors_in_blk))]
        ref_anchors_in_blk$center_point <- rowMeans(ref_anchors_in_blk[,..start_end])

        ref_anchors_in_blk_list <- split(ref_anchors_in_blk, seq(nrow(ref_anchors_in_blk)))

        # Assign a bin to each anchor gene (if possible)
        ref_which_bin <- unlist(lapply(ref_anchors_in_blk_list, function(row){

          # Gene center within a bin
          if(min(abs(ref_center_vec-row$center_point))<(bin_size/2)){
            bin_start <- ref_start_vec[which.min(abs(ref_center_vec-row$center_point))]
            return(bin_start)

            # We assume no gene spans more than 2 bins such that
            # if match the following 2 ifs, these are the only bins possible.

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

        map_dt <- dtw_na(tmp_map_dt = tmp_map_dt, heal_list = heal_list, ref_name = ref_name,
                         alt_name = alt_name, ref_chr = ref_chr, alt_chr = alt_chr, bin_size = bin_size )

        return(map_dt)
      }
      doParallel::stopImplicitCluster()



      }

  }
}
