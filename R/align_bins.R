ref_alt_pair <- function(syn_hits_list,show_which_non_self=TRUE){
  out_list <- lapply(syn_hits_list,function(path){
    filename <- basename(path)
    filename <- gsub(".synHits.txt.gz","",filename)
    which_genomes <- unique(unlist(strsplit(filename,"_vs_")))
    if(show_which_non_self==TRUE){
      return(length(which_genomes)!=1)
    }else{
      return(which_genomes)
    }
  })
  return(out_list)
}

parse_genespace_input <- function(input_dir){

  all_syn_lists <- list.files(input_dir,pattern = "synHits.txt.gz$",recursive = TRUE,full.names = TRUE)
  which_keep <- ref_alt_pair(all_syn_lists,show_which_non_self = TRUE)
  syn_hits_list <- list()
  for(i in 1:length(all_syn_lists)){
    if(which_keep[[i]]==TRUE)  syn_hits_list <- append(syn_hits_list,all_syn_lists[[i]])
  }

  syn_hist_dt <- lapply(syn_hist_list,data.table::fread)
  names(syn_hist_dt) <- ref_alt_pair(syn_hist_list,show_which_non_self = FALSE)

}

align_bins <- function(heal_list,input_dir){

  list.files(input_dir,pattern = "synHits.txt.gz$",recursive = TRUE)
  # Get the map
  syn_blk <- fread(syn_block_path,header = T)
  gnm_oneone <- syn_blk[syn_blk$genome1=="A.halleri" & syn_blk$genome2=="A.lyrata", ]
  one_is_one <- data.frame(gnm_oneone$chr1,gnm_oneone$startBp1,gnm_oneone$endBp1,gnm_oneone$chr2,gnm_oneone$startBp2,gnm_oneone$endBp2,blkID=gnm_oneone$blkID)
  colnames(one_is_one) <- c(paste0(progenitors[1],"_chr"),paste0(progenitors[1],"_start"),paste0(progenitors[1],"_end"),paste0(progenitors[2],"_chr"),paste0(progenitors[2],"_start"),paste0(progenitors[2],"_end"),"blkID")
  map <- one_is_one[order(one_is_one[,1],one_is_one[,2]),]

  # Syntenic anchor genes
  syn_genes <- fread(syn_genes_path,header=T)
  syn_genes <- syn_genes[syn_genes$isAnchor==T,]
  anchors <- data.frame(syn_genes$chr1,syn_genes$start1,syn_genes$end1,syn_genes$chr2,syn_genes$start2,syn_genes$end2,syn_genes$blkID)
  colnames(anchors) <- c(paste0(progenitors[1],"_chr"),paste0(progenitors[1],"_start"),paste0(progenitors[1],"_end"),paste0(progenitors[2],"_chr"),paste0(progenitors[2],"_start"),paste0(progenitors[2],"_end"),"blkID")
  ####
  # SOME BLOCKS ARE MISSING HERE: PROBABLY DIRECTION IS WRONG.
  ####

  #### CHANGE CHROMO NAMES #####
  change2 <- c("chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16")
  names(change2) <- unique(map[,grep("par2_chr",colnames(map))])
  map[,grep("par2_chr",colnames(map))] <- change2[map[,grep("par2_chr",colnames(map))]]
  anchors[,grep("par2_chr",colnames(anchors))] <- change2[anchors[,grep("par2_chr",colnames(anchors))]]



  # Get alignment information
  aln_key <- get_alignment(counts,map,anchors,binSize=binSize,n_cores = n_cores)

}
