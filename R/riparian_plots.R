#' Function to  define a sigmoid style curve between two corrresponding block edges.
#'
#' @param y_vec Points to define along Y. 
#' @param x_start Start point on X.
#' @param x_end End point on X.
#' @param k k param of the sigmoidal function
#'
#' @returns A vector of x axis values to pair with the input y_vec to define the sigmoidal curve.
#'
#' @examples
sigmoid <- function(y, x_start, x_end, k=5) {
  y_min <- min(y)
  y_max <- max(y)
  y0 <- mean(range(y))
  
  x <- x_start + (x_end - x_start) * (1 / (1 + exp(-k * (y - y0))))
  return(x)
}

Dont forget to try and understand why it all works with bottom_x[1], bottom_x[length] but some alt anchors are outside the regions range (on its own chromosome), but this never happens 
with the progenitor. 
plot_riparian <- function(heal_alignment, heal_list, genespace_dir, view_sample = FALSE, output_dir = FALSE, n_threads = 1)
  
  polyploid_samples <- names(heal_alignment)

  progenitors <- gsub("^cn_", "", grep("cn", colnames(heal_alignment[[1]]), value = TRUE))
  
  if(length(progenitors) > 2){
    stop("More than two subgenomes detected. This function only works for two subgenomes at the moment.")
  }
  
  # Get the syntenic regions data table
  path_to_regions_csv <- list.files(path = genespace_dir, 
             pattern = "syntenicRegion_coordinates.csv$", 
             recursive = TRUE, 
             full.names = TRUE)
  regions_dt <- data.table::fread(path_to_regions_csv, quote = "\"")
  # Remove self blocks 
  regions_dt <- regions_dt[!regions_dt$genome1==regions_dt$genome2, ]
  regions_dt <- regions_dt[regions_dt$genome1==progenitors[1]]
  
  # Make a GRanges object
  region_gr_prog1 <- GenomicRanges::GRanges(
    seqnames = regions_dt$chr1,
    ranges = IRanges::IRanges(start = regions_dt$startBp1, end = regions_dt$endBp1),
    region_id = regions_dt$blkID
  )
  
  ## Here I define the dimensions of the plot and the positions of each chromosome
  # First get the approximate length of each chromosome
  subgenome_sizes_list <- foreach::foreach(prog = progenitors)%do%{
    chromosomes <- unique(heal_list[[prog]]$bins$chr)
    size_per_chromo_list <- foreach::foreach(chr = chromosomes)%do%{
      return(max(heal_list[[prog]]$bins$end[heal_list[[prog]]$bins$chr==chr], na.rm = T))
    }
    names(size_per_chromo_list) <- chromosomes
    genome_size <- sum(unlist(size_per_chromo_list))
    output <- list(genome_size=genome_size, size_per_chromo=size_per_chromo_list)
  }
  names(subgenome_sizes_list) <- progenitors
  
  # Make the span a percentage of the longest genome length
  genome_size_subgenomes <- unlist(lapply(subgenome_sizes_list, function(lst){lst$genome_size}))
  longest_genome <- max(genome_size_subgenomes)
  inter_chromosome_space <- 0.05 * longest_genome
  # Get the number of chromosomes and add the interchromosome span times the number of gaps to the genome length
  n_chromo_by_subgenome <- unlist(lapply(subgenome_sizes_list, function(lst){length(lst$size_per_chromo)-1}))
  real_x_spans <- genome_size_subgenomes + n_chromo_by_subgenome * inter_chromosome_space
  names(real_x_spans) <- progenitors
  max_x_span <- max(real_x_spans)
  #####
  ## THERE MIGHT BE SOME USELESS VARIABLES ABOVE
  ######
  # Define y span 
  y_span <- seq(-1, 1, length.out = 300)
  # y_range <- c(0 - 0.5 * y_span, 0 + 0.5 * y_span)
  # y_span <- seq(y_range[1], y_range[2], length.out = 300)
  # 
  offset_list <- foreach::foreach(prog = progenitors)%do%{
    
    chromo_sizes <- unlist(subgenome_sizes_list[[prog]]$size_per_chromo)
    chromo_names <- names(chromo_sizes)
    
    if(length(chromo_sizes)>1){
      
      chromosome_and_spacing <- as.vector(rbind(chromo_sizes, inter_chromosome_space))[-(2 * length(chromo_sizes))]
      names_and_spacing <- as.vector(rbind(chromo_names, "interval"))[-(2 * length(chromo_sizes))]
      
      cum_sum_x_span <- cumsum(chromosome_and_spacing)
      # 0 is start of first chromosome. Then each starts at the end of the preceeding interval. 
      chromo_start_locations <- c(0, cum_sum_x_span[names_and_spacing=="interval"])
      names(chromo_start_locations) <- chromo_names
      
      # We shift the values back to center around 0. Each subgenome is centered around 0 that way. 
      shift_back <- real_x_spans[[prog]]/2
      chromo_start_locations <- chromo_start_locations - shift_back
      
      return(chromo_start_locations)
      
    } else {
      
      start_location <- chromo_sizes / 2
      names(start_location) <- names(chromo_sizes)
      return(start_location)
    }
  }
  names(offset_list) <- progenitors
  
  
  ### Now we can plot each samples
  for(smp in polyploid_samples){
    
    # Create a empty plot.
    p <- ggplot() +
      xlim( -max_x_span/2 + inter_chromosome_space, max_x_span/2 + inter_chromosome_space) +
      ylim(-3, 3) 
    
    # Drawing the chromosomes and their copy number. 
    
    group_vec <- c()
    xmin_vec <-  c()
    xmax_vec <-  c()
    ymin_vec <-  c()
    ymax_vec <-  c()
    
    for(prog in progenitors){
      chr_col_name <- paste0("chr_", prog)
      chromosomes <- unique(heal_alignment[[smp]][[chr_col_name]])
      
      for(chr in chromosomes){
        chromo_dt <- heal_alignment[[smp]][heal_alignment[[smp]][[chr_col_name]]==chr,]
        order_by <- paste0("start_", prog)
        data.table::setorderv(chromo_dt, order_by)
        cn_blocks <- rle(chromo_dt[[paste0("cn_", prog)]])
        block_starts <- c(0, cumsum(cn_blocks$lengths))+1
        
        chromo_dt[[paste0("start_", prog)]][1:block_ends[1]]
        for(b in 1:length(block)){
          group_vec <- c(group_vec, +)
          xmin_vec <-  c(xmin_vec, +)
          xmax_vec <-  c(xmax_vec, +)
          ymin_vec <-  c(ymin_vec, +)
          ymax_vec <-  c(ymax_vec, +)
        }
      }
    }
    
    # Here I define the CN ratio blocks
    blocks <- data.frame(
      group = c("chromosome 1 A.thaliana","chromosome 1 A.thaliana", "chromosome 1 A.arenosa"),
      xmin = c(1, 1, 2),
      xmax = c(4, 4, 6),
      ymin = c(3, 2.4, 7),
      ymax = c(3.5, 2.9, 7.5)
    )
    
    
    
    # Let's define the ribbons (rivers)
    # We work on a alignment ordered by subgenome one 
    ordered_alignment <- heal_alignment[[smp]]
    order_by <- c(paste0("chr_", progenitors[1]), paste0("start_", progenitors[1]))
    data.table::setorderv(ordered_alignment, order_by)
    ### We go through each unique copy number combination AND synteny region. 
    
    ## First we split the alignment in blocks of unique copy number combination
    keep_col <- c(paste0("cn_", progenitors), paste0("chr_", progenitors))
    combination <- do.call(paste, c(ordered_alignment[, ..keep_col], sep = "<sep>"))
    rle_of_blocks <- rle(combination)
    
    start_positions <- c(0, cumsum(rle_of_blocks$lengths))+1
    
    # 
    # ## Then we go through every block and subdivide further by syntenic region
    # # We do the first block outside the loop then added later
    # row_indexes <- 1:end_positions[1]
    # block_1_dt <- ordered_alignment[row_indexes, ]
    # 
    # # Make gene range object
    # block_gr <- GenomicRanges::GRanges(
    #   seqnames = block_1_dt[[paste0("chr_", progenitors[1])]],
    #   ranges = IRanges::IRanges(start = block_1_dt[[paste0("start_",progenitors[1])]], end = block_1_dt[[paste0("end_",progenitors[1])]]),
    #   region_id = block_1_dt[[paste0("id_",progenitors[1])]]
    #   )
    # 
    # # Find overlapping regions
    # regions <- GenomicRanges::findOverlaps(block_gr, region_gr_prog1, select = "first")
    # 
    # # Get the positions of the start and end x values 
    # block_1_ribbon_dts <- foreach::foreach(r = na.omit(unique(regions)))%do%{
    #   
    #   final_block_dt <- block_1_dt[regions == r, ]
    #   top_corners_x <- c(min(final_block_dt[[paste0("start_", progenitors[1])]]), max(final_block_dt[[paste0("end_", progenitors[1])]]))
    #   
    #   # Get order of the second genome region 
    #   direction <- final_block_dt[[paste0("start_", progenitors[2])]][nrow(final_block_dt)] - final_block_dt[[paste0("start_", progenitors[2])]][1] 
    #   if(direction >= 0){
    #     bottom_corners_x <- c(min(final_block_dt[[paste0("start_", progenitors[2])]]), max(final_block_dt[[paste0("end_", progenitors[2])]]))
    #   }else{
    #     bottom_corners_x <- c(max(final_block_dt[[paste0("end_", progenitors[2])]]), min(final_block_dt[[paste0("start_", progenitors[2])]]))
    #   }
    #   
    #   # Now we offset top and bottom positions based on the chromosome. 
    #   prog1_chromo <- final_block_dt[[paste0("chr_", progenitors[1])]][1]
    #   offset_1 <- offset_list[[progenitors[1]]][[prog1_chromo]]
    #   top_corners_x <- top_corners_x + offset_1
    #   
    #   prog2_chromo <- final_block_dt[[paste0("chr_", progenitors[2])]][1]
    #   offset_2 <- offset_list[[progenitors[2]]][[prog2_chromo]]
    #   bottom_corners_x <- bottom_corners_x + offset_2
    #   
    #   # Here I define the edges of the ribbon
    #   ribbon_dt <- data.table::data.table(
    #     y_vec = y_span,
    #     x_start = sigmoid(y_span, bottom_corners_x[1], top_corners_x[1]),
    #     x_end = sigmoid(y_span, bottom_corners_x[2], top_corners_x[2])
    #   )
    #   
    #   return(ribbon_dt)
    # }
    # names(block_1_ribbon_dts) <- regions_dt$blkID[na.omit(unique(regions))]
    # 
    # We go through all other blocks
    ribbon_dt_list <- foreach::foreach(i = 2:length(end_positions))%do%{
  
      row_indexes <- (end_positions[i-1]+1):end_positions[i]
      block_dt <- ordered_alignment[row_indexes, ]
      
      # Make gene range object
      block_gr <- GenomicRanges::GRanges(
        seqnames = block_dt[[paste0("chr_", progenitors[1])]],
        ranges = IRanges::IRanges(start = block_dt[[paste0("start_", progenitors[1])]], end = block_dt[[paste0("end_", progenitors[1])]]),
        region_id = block_dt[[paste0("id_", progenitors[1])]]
      )
# 
#       block_gr_2 <- GenomicRanges::GRanges(
#         seqnames = block_dt[[paste0("chr_", progenitors[2])]],
#         ranges = IRanges::IRanges(start = block_dt[[paste0("start_", progenitors[2])]], end = block_dt[[paste0("end_", progenitors[2])]]),
#         region_id = block_dt[[paste0("id_", progenitors[2])]]
#       )
#       region_gr_prog2 <- GenomicRanges::GRanges(
#         seqnames = regions_dt$chr2,
#         ranges = IRanges::IRanges(start = min(regions_dt$startBp2, regions_dt$endBp2), end = max(regions_dt$startBp2, regions_dt$endBp2)),
#         region_id = regions_dt$blkID
#       )
#       regions_2 <- GenomicRanges::findOverlaps(block_gr_2, region_gr_prog2, select = "first") # If overlapping multiple, only return 1st.

      # Find overlapping regions
      regions <- GenomicRanges::findOverlaps(block_gr, region_gr_prog1, select = "first") # If overlapping multiple, only return 1st. 

      # Get the positions of the start and end x values 
      ribbon_dt_per_block_list <- foreach::foreach(r = na.omit(unique(regions)))%do%{
        
        final_block_dt <- block_dt[regions == r, ]
        top_corners_x <- c(min(final_block_dt[[paste0("start_", progenitors[1])]]), max(final_block_dt[[paste0("end_", progenitors[1])]]))
        
        # Get order of the second genome region 
        direction <- final_block_dt[[paste0("start_", progenitors[2])]][nrow(final_block_dt)] - final_block_dt[[paste0("start_", progenitors[2])]][1] 
        if(direction >= 0){
          #bottom_corners_x <- c(min(final_block_dt[[paste0("start_", progenitors[2])]]), max(final_block_dt[[paste0("end_", progenitors[2])]]))
          bottom_corners_x <- c(final_block_dt[[paste0("start_", progenitors[2])]][1], final_block_dt[[paste0("end_", progenitors[2])]][length(final_block_dt[[paste0("end_", progenitors[2])]])])
        }else if(direction < 0){
          bottom_corners_x <- c(final_block_dt[[paste0("end_", progenitors[2])]][length(final_block_dt[[paste0("end_", progenitors[2])]])], final_block_dt[[paste0("start_", progenitors[2])]][1])
          #bottom_corners_x <- c(max(final_block_dt[[paste0("end_", progenitors[2])]]), min(final_block_dt[[paste0("start_", progenitors[2])]]))
        }
        
        # Now we offset top and bottom positions based on the chromosome. 
        prog1_chromo <- final_block_dt[[paste0("chr_", progenitors[1])]][1]
        offset_1 <- offset_list[[progenitors[1]]][[prog1_chromo]]
        top_corners_x <- top_corners_x + offset_1
        
        prog2_chromo <- final_block_dt[[paste0("chr_", progenitors[2])]][1]
        offset_2 <- offset_list[[progenitors[2]]][[prog2_chromo]]
        bottom_corners_x <- bottom_corners_x + offset_2
        
        # Here I define the edges of the ribbon
        ribbon_dt <- data.table::data.table(
          y_vec = y_span,
          x_start = sigmoid(y_span, bottom_corners_x[1], top_corners_x[1]),
          x_end = sigmoid(y_span, bottom_corners_x[2], top_corners_x[2])
        )
        
        return(ribbon_dt)
      }
      names(ribbon_dt_per_block_list) <- regions_dt$blkID[na.omit(unique(regions))]
      return(ribbon_dt_per_block_list)
    }
    
    ribbon_dt_list[[length(ribbon_dt_list)+1]] <- block_1_ribbon_dts
    ribbon_dt_list <- do.call(c, ribbon_dt_list) # We combine all lists into one big list. Each entry is a ribbon defining data table.
    
    ### Add the ribbons to the plot
    for (i in seq_along(ribbon_dt_list)) {
      p <- p + ggplot2::geom_ribbon(
          data = as.data.frame(ribbon_dt_list[[i]]),
          aes(y = y_vec, xmin = x_start, xmax = x_end),
          fill = "grey50", alpha = 0.5
      )
    }
    
    
    
    # Here I define the CN ratio blocks
    blocks <- data.frame(
      group = c("chromosome 1 A.thaliana","chromosome 1 A.thaliana", "chromosome 1 A.arenosa"),
      xmin = c(1, 1, 2),
      xmax = c(4, 4, 6),
      ymin = c(3, 2.4, 7),
      ymax = c(3.5, 2.9, 7.5)
    )
    
    
    
    # Here I define the edges of the ribbon
    ribbon_df <- data.frame(
      y = y_vals,
      x_left = sigmoid(y_vals, 1, 2),
      x_right = sigmoid(y_vals, 4, 6)
    )
    
    
    # Plot
    ggplot2::ggplot() +
      ggplot2::geom_rect(data = blocks, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = group), color = "black", alpha = 0.5) +
      ggplot2::geom_ribbon(data = ribbon_df, aes(y = y, xmin = x_left, xmax = x_right),
                  fill = "grey50", alpha = 0.5) +
      ggplot2::theme_minimal()
  }