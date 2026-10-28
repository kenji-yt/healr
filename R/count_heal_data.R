#' Parse feature count output for bins and merge to gc and mappability.
#'
#' @param feature_count_list Feature counts output list.
#' @param gc_dt GC content data table.
#' @param map_dt Mappability  data table.
#' @param sample_names Sample names.
#'
#' @return A data table with GC, mappability and count per bin.
#'
parse_counts_bins <- function(feature_count_list, gc_dt, map_dt, sample_names) {
  counts_dt <- data.table::as.data.table(feature_count_list$counts)
  colnames(counts_dt) <- sample_names

  gene_ID <- feature_count_list$annotation$GeneID
  chr <- feature_count_list$annotation$Chr
  start <- feature_count_list$annotation$Start
  end <- feature_count_list$annotation$End
  length <- feature_count_list$annotation$Length
  
  bin_size <- stats::median(feature_count_list$annotation$Length)

  anno_dt <- data.table::data.table(chr = chr[length == bin_size], start = start[length == bin_size] - 1, end = end[length == bin_size])

  out_dt <- cbind(anno_dt, counts_dt[length == bin_size, ])
  data.table::setkey(out_dt, chr, start)
  data.table::setorder(out_dt, chr, start)

  out_dt <- merge(gc_dt, out_dt, by = c("chr", "start"))
  out_dt <- merge(map_dt, out_dt, by = c("chr", "start"))
  data.table::setkey(out_dt, chr, start, gc_content, mappability)

  return(out_dt)
}


#' Counts reads and combines results with GC and mappability information.
#'
#' @param input_dir A healr input directory (see details).
#' @param n_threads Number of threads to use with featureCounts ('1' by default)
#' @param paired_end Logical: Is the data paired end.
#' @param full_output Logical: Do you want to also get the full featureCounts output ('FALSE' by default). This can be usefull for debugging. 
#' @param unknown_reads Logical: Also count reads that could not be assigned to a subgenome ('TRUE' by default). Their BAM file must be in a subdirectory 'unknown' of each sample's subgenome directory and finish with ".unk.bam". Adds one column per sample to the ‘bin_info’ table with the number of unknown reads in each bin (column name '<sample>_unknown').
#'
#' @return A heal list object i.e. a list with one element per progenitor containing a 'bins' data table with counts in bins for each sample and GC and mappability for each bin, and ‘bin_info’ data table with GC content and mappability (and unknown reads optionally) for all bins, which is kept when bins are filtered.
#' @export
#'
#' @importFrom foreach %do%
#' @importFrom utils write.table
count_heal_data <- function(input_dir, n_threads = 1, paired_end, full_output = FALSE, unknown_reads = TRUE) {
  
  prog_dir <- list.files(path = input_dir, pattern = "progenitors", full.names = TRUE)
  if(length(prog_dir)==0){
    stop("No directory named progenitors in the input directory.")
  }else if(length(prog_dir)>1){
    stop("More than one directory named progenitors in the input directory.")
  }
  
  poly_dir <- list.files(path = input_dir, pattern = "polyploids", full.names = TRUE)
  if(length(poly_dir)==0){
    stop("No directory named polyploids in the input directory.")
  }else if(length(poly_dir)>1){
    stop("More than one directory named polyploids in the input directory.")
  }
  
  progenitors <- list.files(path = prog_dir)

  # Check that every sample has its BAM file(s) for every subgenome before counting anything
  sample_dirs <- list.dirs(path = poly_dir, full.names = TRUE, recursive = FALSE)
  for (prog in progenitors) {
    bam_names <- basename(dirname(dirname(list.files(file.path(sample_dirs, prog), pattern = "\\.bam$", full.names = TRUE))))
    if (anyDuplicated(bam_names) || length(bam_names) != length(sample_dirs)) {
      stop(paste0("Each sample directory in ", poly_dir, " must contain exactly one .bam file in a subdirectory named ", prog, "."))
    }
    if (unknown_reads == TRUE) {
      unk_names <- basename(dirname(dirname(dirname(list.files(file.path(sample_dirs, prog, "unknown"), pattern = "\\.unk\\.bam$", full.names = TRUE)))))
      if (anyDuplicated(unk_names) || length(unk_names) != length(sample_dirs)) {
        stop(paste0("Each sample directory in ", poly_dir, " must contain exactly one '.unk.bam' file in ", file.path(prog, "unknown"), ".\n",
                    "If you don't have the 'unk.bam' files from EAGLE-RC or do not wish to consider unknown reads, pass 'unknown_reads = FALSE' as argument."))
      }
    }
  }
  
  counts <- foreach::foreach(prog = progenitors) %do% {
    # get paths ----------------

    cat(paste0("Counting reads in bins for progenitor or subgenome ", prog, ".", "\n"))

    current_prog_dir <- file.path(prog_dir, prog)
    
    gc_path <- list.files(path = current_prog_dir, pattern = "gc.bed$", full.names = TRUE)
    map_path <- list.files(path = current_prog_dir, pattern = "mappability.bed$", full.names = TRUE)
    bins_path <- list.files(path = current_prog_dir, pattern = "bins.bed$", full.names = TRUE)

    if (length(gc_path) != 1 || length(map_path) != 1 || length(bins_path) != 1) {
      stop(paste0("Not exactly one mappability, GC content or bins file for progenitor ", prog, ". Exiting!"))
    }


    # load data --------------------

    gc_dt <- data.table::fread(gc_path, select = c(1, 2, 5))
    colnames(gc_dt) <- c("chr", "start", "gc_content")
    data.table::setkey(gc_dt, chr, start, gc_content)
    data.table::setorder(gc_dt, chr, start)

    map_dt <- data.table::fread(map_path, select = c(1, 2, 4))
    colnames(map_dt) <- c("chr", "start", "mappability")
    data.table::setkey(map_dt, chr, start, mappability)
    data.table::setorder(map_dt, chr, start)

    bins_dt <- data.table::fread(bins_path, select = c(1, 2, 3))
    colnames(bins_dt) <- c("chr", "start", "end")
    data.table::setkey(bins_dt, chr, start)
    data.table::setorder(bins_dt, chr, start)


    # make SAF -------------------

    bins_saf <- data.frame(
      GeneID = paste0(bins_dt$chr, "_", bins_dt$start + 1),
      Chr = bins_dt$chr,
      Start = bins_dt$start + 1,
      End = bins_dt$end,
      Strand = "+"
    )

    anno_bins <- paste0(dirname(bins_path), "/bins.saf")

    withr::with_options(list(scipen = 999), write.table(bins_saf, file = anno_bins, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE))


    # count in bins --------------

    bam_paths <- list.files(file.path(sample_dirs, prog), pattern = "\\.bam$", full.names = TRUE)
    sample_names <- basename(dirname(dirname(bam_paths)))
    feature_count_list <- Rsubread::featureCounts(bam_paths, annot.ext = anno_bins, isPairedEnd = paired_end, nthreads = n_threads, allowMultiOverlap = TRUE)
    count_dt <- parse_counts_bins(feature_count_list = feature_count_list, gc_dt = gc_dt, map_dt = map_dt, sample_names = sample_names)
    
    out <- list(bins = count_dt, bin_info = count_dt[, c("chr", "start", "end", "gc_content", "mappability"), with = FALSE])
    
    # count unknown reads -----------
    if (unknown_reads == TRUE) {
      unk_paths <- list.files(file.path(sample_dirs, prog, "unknown"), pattern = "\\.unk\\.bam$", full.names = TRUE)
      unk_names <- basename(dirname(dirname(dirname(unk_paths))))
      unk_count_list <- Rsubread::featureCounts(unk_paths, annot.ext = anno_bins, isPairedEnd = paired_end, nthreads = n_threads, allowMultiOverlap = TRUE)
      unk_dt <- parse_counts_bins(feature_count_list = unk_count_list, gc_dt = gc_dt, map_dt = map_dt, sample_names = unk_names)
      unk_dt <- unk_dt[, c("chr", "start", unk_names), with = FALSE]
      data.table::setnames(unk_dt, unk_names, paste0(unk_names, "_unknown"))
      out$bin_info <- merge(out$bin_info, unk_dt, by = c("chr", "start"))
    }
    
    if (full_output == TRUE) {
      out$feature_counts_output <- feature_count_list
    }
    return(out)
    
  }

  names(counts) <- progenitors

  return(counts)
}


utils::globalVariables(c("prog", "chr", "start", "end", "GeneID", "gc_content", "mappability"))
