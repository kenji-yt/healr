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
sigmoid <- function(y, x_start, x_end, k=7) {
  y_min <- min(y)
  y_max <- max(y)
  y0 <- mean(range(y))
  
  x <- x_start + (x_end - x_start) * (1 / (1 + exp(-k * (y - y0))))
  return(x)
}

#Dont forget to try and understand why it all works with bottom_x[1], bottom_x[length] but some alt anchors are outside the regions range (on its own chromosome), but this never happens 
#with the progenitor. 
### PARALELIZE!!!
#' Plot riparian style plots with copy number
#'
#' @param alignment A heal alignment object created with get_heal_alignment().
#' @param heal_list List in heal format (such as output from count_heal_data()).
#' @param genespace_dir Path to a directory containing the syntenicRegion_coordinates.csv output file from GENESPACE.
#' @param output_dir The name of a directory to write all plots to. Will create one if nonexistent.
#' @param n_threads Number of threads to use ('1' by default).
#' @param colour_scales A list with one entry per subgenome. Each entry must be named after a subgenome and contain a vector with the starting and ending colour value of the range. Defaults to red and green colour ranges.
#' @param theme Background settings. Options are 'dark' or 'light'. Default value is 'light'. 
#' 
#' @returns Nothing. Creates riparian style plots with copy number information. 
#' @export
#'
#' @examples
plot_riparian <- function(heal_alignment, heal_list, genespace_dir, output_dir = FALSE, n_threads = 1, colour_scales = FALSE, theme = 'light'){
  
  if(theme != "light" && theme != "dark"){
    stop("Invalid input for 'theme'. Please input 'light' or 'dark'. Exiting..")
  }
  polyploid_samples <- names(heal_alignment)

  progenitors <- gsub("^cn_", "", grep("cn", colnames(heal_alignment[[1]]), value = TRUE))
  
  if(length(progenitors) > 2){
    stop("More than two subgenomes detected. This function only works for two subgenomes at the moment. Exiting..")
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
  
  # Make the inter chromosome gap a percentage of the longest genome length
  genome_size_subgenomes <- unlist(lapply(subgenome_sizes_list, function(lst){lst$genome_size}))
  longest_genome <- max(genome_size_subgenomes)
  inter_chromosome_space <- 0.05 * longest_genome
  # Get the number of chromosomes and add the interchromosome gap times the number of gaps to get the genome length(s)
  n_chromo_by_subgenome <- unlist(lapply(subgenome_sizes_list, function(lst){length(lst$size_per_chromo)}))
  real_x_spans <- genome_size_subgenomes + (n_chromo_by_subgenome-1) * inter_chromosome_space
  names(real_x_spans) <- progenitors
  max_x_span <- max(real_x_spans)
  
  # Define y span 
  y_span <- seq(-1, 1, length.out = 300)
 
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
  
  # Y location of each subgenome
  direction <- c(1, -1) 
  names(direction) <- progenitors
  # Edges (absolution from 0)
  cn_block_edges <- list(c(1, 1.2), c(1.3, 1.5), c(1.7, 1.9), c(2.1, 2.3))
  
  ### Now we can plot each samples
  for(smp in polyploid_samples){
    
    # Create an empty plot.
    plot <- ggplot2::ggplot() +
      ggplot2::xlim(-max_x_span/2 - inter_chromosome_space, max_x_span/2 + inter_chromosome_space) +
      ggplot2::ylim(-2.6, 2.6) +
      ggplot2::theme_void() +
      ggplot2::labs(title = smp) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold")
      )
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
        
        if(length(block_starts)>1){
          
          for(i in 1:(length(block_starts)-1)){
            
            start_vec <- chromo_dt[[paste0("start_", prog)]][block_starts[i]:(block_starts[i+1]-1)]
            end_vec <- chromo_dt[[paste0("end_", prog)]][block_starts[i]:(block_starts[i+1]-1)]
            cn <- unique(chromo_dt[[paste0("cn_", prog)]][block_starts[i]:(block_starts[i+1]-1)])
            
            if(!is.na(cn)){
              
              if(cn>0){
                
                offset <- offset_list[[prog]][[chr]]
                start_loc <- min(start_vec) + offset
                end_loc <- max(end_vec) + offset
                
                if(cn > 4){
                  bloc_center <- mean(start_loc, end_loc)
                  y_position <- 2.5 * direction[[prog]]
                  plot <- plot + ggplot2::annotate("text", x = bloc_center, y = y_position, label = "+", size = 6) 
                  cn <- 4
                }
                
                for(i in 1:cn){
                  
                  group_vec <- c(group_vec, paste(chr, prog))
                  xmin_vec <-  c(xmin_vec, start_loc)
                  xmax_vec <-  c(xmax_vec, end_loc)
                  
                  ymin <- cn_block_edges[[i]][1] * direction[[prog]]
                  ymax <- cn_block_edges[[i]][2] * direction[[prog]]
                  
                  ymin_vec <-  c(ymin_vec, ymin)
                  ymax_vec <-  c(ymax_vec, ymax)
                }
              }
            }
          }
        }
      }
    }
    
    
    # Here I define the CN ratio blocks
    blocks_dt <- data.frame(
      Chromosome = group_vec,
      xmin = xmin_vec,
      xmax = xmax_vec,
      ymin = ymin_vec,
      ymax = ymax_vec
    )
    
    blocks_dt$Chromosome <- factor(
      blocks_dt$Chromosome,
      levels = c(sort(unique(grep(progenitors[1], blocks_dt$Chromosome, value = TRUE))), sort(unique(grep(progenitors[2], blocks_dt$Chromosome, value = TRUE))))
    )

    if(is.logical(colour_scales)){
      if(colour_scales==TRUE){
        stop("Invalid entry for 'colour_scales'. Please input a list with one entry per subgenome. Each entry must be named after a subgenome and contain a vector with the starting and ending colour value of the range.")
      }
      subgenome_1_scale <- c("lightgreen", "darkgreen")
      subgenome_2_scale <- c("lightcoral", "darkred")
      colour_scales <- list(subgenome_1_scale, subgenome_2_scale)
      names(colour_scales) <- names(n_chromo_by_subgenome) 
    }else{
      if(length(intersect(names(colour_scales), names(n_chromo_by_subgenome))) == length(progenitors)){
        if(sum(names(colour_scales) == names(n_chromo_by_subgenome)) != length(progenitors)){
          colour_scales <- list(colour_scales[[2]], colour_scales[[1]])
          names(colour_scales) <- names(n_chromo_by_subgenome)
        }
      }else{
        stop("Invalid entry for 'colour_scales'. Please input a list with one entry per subgenome. Each entry must be named after a subgenome and contain a vector with the starting and ending colour value of the range.")
      }
    }
    
    plot <- plot + ggplot2::geom_rect(data = blocks_dt, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Chromosome), color = "#272727") +
      ggplot2::scale_fill_manual(values = c(
        scales::gradient_n_pal(colour_scales[[progenitors[1]]])(seq(0, 1, length.out = n_chromo_by_subgenome[[progenitors[1]]])),
        scales::gradient_n_pal(colour_scales[[progenitors[2]]])(seq(0, 1, length.out = n_chromo_by_subgenome[[progenitors[2]]])))
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
    
    # We go through all blocks
    ribbon_dt_list <- foreach::foreach(i = 1:(length(start_positions)-1))%do%{
  
      row_indexes <- start_positions[i]:(start_positions[i+1]-1)
      block_dt <- ordered_alignment[row_indexes, ]
      
      # Make gene range object
      block_gr <- GenomicRanges::GRanges(
        seqnames = block_dt[[paste0("chr_", progenitors[1])]],
        ranges = IRanges::IRanges(start = block_dt[[paste0("start_", progenitors[1])]], end = block_dt[[paste0("end_", progenitors[1])]]),
        region_id = block_dt[[paste0("id_", progenitors[1])]]
      )

      # Find overlapping regions
      regions <- GenomicRanges::findOverlaps(block_gr, region_gr_prog1, select = "first") # If overlapping multiple, only return 1st. 

      # Get the positions of the start and end x values 
      ribbon_dt_per_block_list <- foreach::foreach(r = na.omit(unique(regions)))%do%{
        
        final_block_dt <- block_dt[regions == r, ]
        top_corners_x <- c(min(final_block_dt[[paste0("start_", progenitors[1])]]), max(final_block_dt[[paste0("end_", progenitors[1])]]))
        
        # Get order of the second genome region 
        orientation <- final_block_dt[[paste0("start_", progenitors[2])]][nrow(final_block_dt)] - final_block_dt[[paste0("start_", progenitors[2])]][1] 
        if(orientation >= 0){
          bottom_corners_x <- c(final_block_dt[[paste0("start_", progenitors[2])]][1], final_block_dt[[paste0("end_", progenitors[2])]][length(final_block_dt[[paste0("end_", progenitors[2])]])])
        }else if(orientation < 0){
          bottom_corners_x <- c(final_block_dt[[paste0("end_", progenitors[2])]][length(final_block_dt[[paste0("end_", progenitors[2])]])], final_block_dt[[paste0("start_", progenitors[2])]][1])
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
    
    ribbon_dt_list <- do.call(c, ribbon_dt_list) # We combine all lists into one big list. Each entry is a ribbon defining data table.
    
    ### Add the ribbons to the plot
    for (i in seq_along(ribbon_dt_list)) {
      plot <- plot + ggplot2::geom_ribbon(
          data = as.data.frame(ribbon_dt_list[[i]]),
          ggplot2::aes(y = y_vec, xmin = x_start, xmax = x_end),
          fill = "grey50", alpha = 0.5
      )
    }
    
    if(theme=="dark"){
      plot <- plot +
        ggplot2::theme(
          plot.background = ggplot2::element_rect(fill = "black"),
          panel.background = ggplot2::element_rect(fill = "black"),
          #panel.grid.major = ggplot2::element_line(color = "gray30"),
          #panel.grid.minor = ggplot2::element_line(color = "gray20"),
          text = ggplot2::element_text(color = "white"),
          #axis.text = ggplot2::element_text(color = "white"),
          #axis.title = ggplot2::element_text(color = "white"),
          #plot.title = ggplot2::element_text(color = "white", size = 16, face = "bold"),
          legend.background = ggplot2::element_rect(fill = "black"),
          legend.key = ggplot2::element_rect(fill = "black"),
          legend.text = ggplot2::element_text(color = "white"),
          legend.title = ggplot2::element_text(color = "white")
        )
    }
    
    print(plot)
    
    if(!is.logical(output_dir)){
      dir.create(file.path(output_dir))
      ggplot2::ggsave(file=paste0(output_dir, "/", smp, "_riparian.svg"), plot=plot, width=20, height=8)
      ggplot2::ggsave(file=paste0(output_dir, "/", smp, "_riparian.pdf"), plot=plot, width=20, height=8)
    }
  }
}