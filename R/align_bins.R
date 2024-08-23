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


  return(list(syn_hits=syn_hits_dt,block_coord=blk_coord_dt))
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


      blks <- map_dt_list$block_coord[map_dt_list$block_coord$genome1==ref_gnm & map_dt_list$block_coord$genome2==alt_gnm,]

      doParallel::registerDoParallel(n_cores)

      foreach::foreach(i=(1:nrow(blks)))%do%{

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

        ref_start_vec <- seq(ref_bin_start,ref_bin_end,bin_size*sign(ref_bin_end-ref_bin_start))
        ref_chromo_vec <- rep(ref_chr,length(ref_start_vec))

        ref_center_vec <- ref_start_vec+(bin_size/2)



        anchors_in_blk <- anchors_dt[anchors_dt$blkID==blk_id,]

        ref_anchor_centers <- rowMeans(anchors_in_blk[,c("start1","end1")])

        which_bin_ref <- ref_start_vec[sapply(ref_anchor_centers,function(x){which.min(abs(ref_center_vec-x))})]

        alt_anchor_centers <- rowMeans(anchors_in_blk[,c("start2","end2")])
        names(alt_anchor_centers) <- which_bin_ref

        map_list <- sapply(ref_start_vec,function(bin){list(ref_start=bin,alt_center=mean(alt_anchor_centers[as.character(bin)]))}, simplify = FALSE)
        map_df <- do.call(rbind, lapply(map_list,as.data.frame))

        # Infer position of bins with no anchor genes.
        map_df$alt_center <- replace_na(map_df$alt_center,alt_start = alt_start,alt_end = alt_end)

        alt_chromo_vec <- rep(alt_chr,length(ref_start_vec))


        alignment <- data.frame(chr=ref_chromo_vec,start=map_df$ref_start,map_chr=alt_chromo_vec,map_pos=map_df$alt_center)
        return(alignment)


      }
      doParallel::stopImplicitCluster()
      }

  }
}
