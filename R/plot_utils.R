#' Plot counts and copy number for multiple samples together
#'
#' @param heal_list List in heal format (such as output from count_heal_data()).
#' @param view_samples A vector of sample names to be plotted. Samples will be plotted in the order provided. Default value is "all", which shows all samples. 
#' @param output_dir The name of a directory to write all plots to. Will create one if nonexistent.
#' @param prog_ploidy of the progenitors (Assumed to be equal. '2' by default).
#' @param plot_cn Logical: plot a line indicating infered copy number ('FALSE' by default in CN has been estimated).
#' @param add_bins Logical: plot counts for each bin ('TRUE' by default; normalized in plot_cn=TRUE).
#' @param color_map A vector of colors for each progenitor. If "FALSE" the colors are choosen using viridis(). 
#' @param method Which method was used to assign a copy number to each segment (if plot_cn=TRUE and add_bins!=FALSE). It should be the same method as the one used in get_copy_number() ('global', 'local' or 'manual'. 'global' by default). 
#' @param average The method to compute average to normalize counts (if plot_cn=TRUE and add_bins!=FALSE). It should be the same method as the one used in get_copy_number() ('median' or 'mean'. 'median' by default). 
#' @param average_list Same as in get_copy_number(). 
#' @param return_plot Logical: return the ggplot2 object ('FALSE' by default).
#' @param chr_label_size Size of the chromosome labels as input to ggplot2::geom_text size argument ('4' by default).
#' @param cn_label_size Size of the copy number labels as input to ggplot2::geom_text size argument ('3' by default).
#' @param size_x_axis_title Size of the x axis title as input to ggplot2::element_text size argument ('14' by default).
#' @param size_y_axis_title Size of the x axis title as input to ggplot2::element_text size argument ('14' by default).
#' @param sample_name_size Size of the sample names as input to ggplot2::geom_text size argument ('3' by default).
#' @param bin_point_alpha The transparency of the normalized bin count values as input to ggplot2::goem_point alpha argument ('0.4' by default). 
#' @param bin_point_size Size of the normalized bin count values as input to ggplot2::goem_point size argument ('2' by default).
#' @param cn_line_width Thickness of the line showing the copy number as input to ggplot2::geom_path linewidth argument ('2' by default).
#' @param device Plot device (as argument to ggplot2::ggsave) ('pdf' by default).
#' @param ... Any arguments you wish to pass to ggplot2::ggsave().
#' 
#' @return Either nothing of a ggplot object showing several samples together. The plot is always printed. 
#' @export
#' 
#' @examples
plot_all_bins <- function(heal_list, view_samples = "all", output_dir = FALSE,
                          prog_ploidy = 2, plot_cn = FALSE, add_bins = TRUE, color_map = FALSE,
                          method = "global", average = "median", average_list = FALSE, return_plot = FALSE,
                          chr_label_size = 4, cn_label_size = 3, size_x_axis_title = 14, 
                          size_y_axis_title = 14, sample_name_size = 3,  bin_point_alpha = 0.4,
                          bin_point_size = 2, cn_line_width = 2, device = "pdf", ...){
  
  # Check input to see if it's appropriate
  if(!is.logical(add_bins)){
    stop("add bins must be either TRUE or FALSE (logical). Exiting..")
  }
  cn_is_null <- is.null(unlist(lapply(heal_list, function(list) {
    list$CN
  })))
  if (cn_is_null && plot_cn == TRUE) {
    cat("ERROR: plot_cn set to 'TRUE' but no CN data table found in heal_list. Setting plot_cn to 'FALSE'. \n")
    plot_cn <- FALSE
  }
  
  # Get the averages for normalization
  if(method == "global"){
    average_list <- get_sample_stats(heal_list, method = average)
    sample_names <- names(average_list)
    
  }else if(method == "local"){
    average_list <- get_sample_stats(heal_list, method = paste0("local_", average))
    sample_names <- names(average_list)
    
  }else if (method == "manual" & average_list != FALSE){
    sample_names <- names(average_list)
  }
  
  
  progenitors <- names(heal_list)
  
  # Get the colour map
  if (isFALSE(color_map)) {
    color_map <- viridis::viridis(length(progenitors))
  } else if (length(color_map) != length(progenitors)) {
    stop("Custom color_map is not of correct length. It should match the number of subgenomes. Exiting..")
  }
  names(color_map) <- progenitors

  # Define which samples the users wishes on the plot
  if(length(view_samples)==1){
    if(view_samples=="all"){
      view_samples <- sample_names 
    }
  }
  
  if (length(intersect(sample_names, view_samples)) == 0) {
    stop("Sample names not recognized for viewing of counts or coverage. Exiting..")
  }
  
  samples <- rev(view_samples)
       
  cat(paste0("Plotting copy number for ", paste(samples, collapse = ", "), ".\n"))
  
  if (output_dir != FALSE) {
    if(dir.exists(output_dir)){
      cat(paste0("Saving plot to ", output_dir, ".", "\n"))
    }else{
      stop(paste0("Output directory ", output_dir, " does not exist. Exiting.."))
    }
  }

    
  
  # Get y position of each sample
  sample_y_starts <- 0:(length(samples)-1)*8
  names(sample_y_starts) <- samples
  
  
  ## Here we define the dimensions of the plot and the positions of each chromosome
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
  n_chromo_by_subgenome <- unlist(lapply(subgenome_sizes_list, function(lst){length(lst$size_per_chromo)-1}))
  real_x_spans_per_subG <- genome_size_subgenomes + n_chromo_by_subgenome * inter_chromosome_space
  names(real_x_spans_per_subG) <- progenitors
  total_x_span <- sum(real_x_spans_per_subG) + (inter_chromosome_space * 2) * (length(progenitors)-1) # add two interchromosome gaps between subgenomes
  
  # Create an empty plot.
  plot <- ggplot2::ggplot() +
    ggplot2::xlim(-inter_chromosome_space, total_x_span+inter_chromosome_space*0.1) +
    ggplot2::ylim(-2, max(sample_y_starts)+8) +
    ggplot2::theme_void()
  
  # Add sample names
  smp_names_df <- data.frame(x = 0, y = sample_y_starts + prog_ploidy*length(progenitors) + prog_ploidy, label = names(sample_y_starts))
  plot <- plot + 
    ggplot2::geom_text(data = smp_names_df,
                       ggplot2::aes(x = x, y = y, label = label, fontface = "bold"),
                       hjust = 0,           # align text to the right
                       size = sample_name_size)
  
  # Add y axes
  y_ends <- sample_y_starts + prog_ploidy*length(progenitors) 
  y_axes_df <- data.frame(x_start = 0, x_end = 0, y_start = sample_y_starts, y_end = y_ends)
  plot <- plot + 
    ggplot2::geom_segment(data = y_axes_df,
                          ggplot2::aes(x = x_start, xend = x_end, y = y_start, yend = y_end))
  
  # Add ticks 
  tick_location <- sample_y_starts
  for(i in 1:(prog_ploidy*length(progenitors))){
    tick_location <- c(tick_location, sample_y_starts + i)
  }
  y_axes_ticks_df <- data.frame(x_start = -inter_chromosome_space*0.6, x_end = 0, y_start = tick_location, y_end = tick_location)
  
  tick_labels <- rep(0:(prog_ploidy*length(progenitors)), length(sample_y_starts))
  y_ticks_labels_df <- data.frame(x = -inter_chromosome_space*0.8, y = sort(tick_location), label = tick_labels)
  
  plot <- plot + 
    ggplot2::geom_segment(data = y_axes_ticks_df,
                          ggplot2::aes(x = x_start, xend = x_end, y = y_start, yend = y_end))+
    ggplot2::geom_text(data = y_ticks_labels_df,
                       ggplot2::aes(x = x, y = y, label = label),
                       hjust = 1, 
                       size = cn_label_size)
  
  # Define progenitor specific offsets
  prog_offsets <- c(0, real_x_spans_per_subG + (inter_chromosome_space * 2))
  prog_offsets <- prog_offsets[1:length(progenitors)]
  names(prog_offsets) <- progenitors
  
  # Get the starting position of each chromosome.
  x_start_vec_list <- foreach::foreach(prog = progenitors)%do%{
    chr_vec <- unlist(subgenome_sizes_list[[prog]]$size_per_chromo)
    c(0, cumsum(chr_vec + inter_chromosome_space) + 1)
    start_pos <- c(0, cumsum(chr_vec + inter_chromosome_space) + 1)
    start_pos <- start_pos[1:length(chr_vec)]
    start_pos <- start_pos + prog_offsets[[prog]]
    
    names(start_pos) <- names(chr_vec)
    return(start_pos)
  }
  names(x_start_vec_list) <- progenitors
  
  # Add chromosome ranges, ticks and names 
  # Ranges
  x_starts <- unlist(x_start_vec_list)
  chromo_sizes <- unlist(lapply(subgenome_sizes_list, function(list){list$size_per_chromo}))
  x_ends <- x_starts + chromo_sizes
  chromo_axes_df <- data.frame(x_start = x_starts, x_end = x_ends, y_start = -1, y_end = -1)
  # Ticks
  tick_x_location <- unlist(apply(chromo_axes_df, 1, function(row){
    seq(row[1], row[2], by=5000000)
  }))
  chromo_ticks_df <- data.frame(x_start = tick_x_location, x_end = tick_x_location, y_start = -1, y_end = -1.4)
  # Names
  label_positions <- rowMeans(chromo_axes_df[, c("x_start", "x_end")])
  chromo_names <- unlist(lapply(subgenome_sizes_list, function(list){names(list$size_per_chromo)}))
  chromo_labels_df <- data.frame(x = label_positions, y = -2, label = chromo_names)
  
  plot <- plot + 
    ggplot2::geom_segment(data = chromo_axes_df,
                          ggplot2::aes(x = x_start, xend = x_end, y = y_start, yend = y_end))+ 
    ggplot2::geom_segment(data = chromo_ticks_df,
                          ggplot2::aes(x = x_start, xend = x_end, y = y_start, yend = y_end),
                          linewidth = 0.1)+ 
    ggplot2::geom_text(data = chromo_labels_df,
                       ggplot2::aes(x = x, y = y, label = label),
                       hjust = 0.5, 
                       size = chr_label_size)
  
  # Add dotted lines for CN
  for(i in 1:nrow(chromo_axes_df)){
    x_start <- chromo_axes_df$x_start[i]
    x_end <- chromo_axes_df$x_end[i]
    y_starts <- y_axes_ticks_df$y_start
    y_ends <- y_axes_ticks_df$y_end
    
    dotted_df <- data.frame(x_start = x_start, x_end = x_end, y_start = y_starts, y_end = y_ends)
    plot <- plot + 
      ggplot2::geom_segment(data = dotted_df,
                            ggplot2::aes(x = x_start, xend = x_end, y = y_start, yend = y_end),
                            linewidth = 0.2,            
                            alpha = 0.8,        
                            linetype = "dotted" 
      )
  }
  
  # Add a legend and x and y labels 
  plot <- plot + ggplot2::geom_point(data = data.frame(name = factor(names(color_map), levels = names(color_map))), 
                                     ggplot2::aes(x = Inf, y = Inf, color = name), 
                                     alpha = 2) + 
    ggplot2::scale_color_manual(values = color_map, 
                                labels = names(color_map), 
                                name = "Subgenomes",
                                guide = ggplot2::guide_legend(override.aes = list(size = 5))  # make legend points bigger
    ) +
    ggplot2::theme(
      legend.position = "right",
      legend.box.margin = ggplot2::margin(0, 0, 0, -10)  # negative right margin pulls legend closer
    ) +
    ggplot2::labs(
      x = "Chromosomes",
      y = "Copy Number"
    ) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(color = "black", size = size_x_axis_title, 
                                           vjust = 4),
      axis.title.y = ggplot2::element_text(color = "black", size = size_y_axis_title, angle = 90,
                                           vjust = -4)
    )
  
  
  # Add bins and copy numbers to the plot
  for(prog in progenitors){
    
    current_chromo <- unique(heal_list[[prog]]$bins$chr)
    
    for(smp in samples){
      
      # Set a max allowed CN
      max_value <- prog_ploidy*(length(progenitors)+1)
      
      # Normalize and filter
      normalized <- heal_list[[prog]]$bins[[smp]] / average_list[[smp]][[prog]] * prog_ploidy
      normalized[normalized > max_value] <- NA
      
      cn <- heal_list[[prog]]$CN[[smp]]
      which_plus <- cn > max_value
      cn[which_plus] <- max_value+0.1
      
      for(chr in current_chromo){
        
        which_rows <- heal_list[[prog]]$bins$chr == chr
        x_vec <- heal_list[[prog]]$bins$start[which_rows] + x_start_vec_list[[prog]][[as.character(chr)]]
        y_pts <- normalized[which_rows] + sample_y_starts[[smp]]
        
        point_df <- data.frame(x=x_vec, y=y_pts)
        
        if(plot_cn==TRUE){
          y_line <- cn[which_rows] + sample_y_starts[[smp]]
          line_df <- data.frame(x=x_vec, y=y_line)
          
          if(add_bins==TRUE){
            
            plot <- plot + 
              ggplot2::geom_point(data = point_df, ggplot2::aes(x = x, y = y), size = bin_point_size, color = "black", alpha=bin_point_alpha) +
              ggplot2::geom_path(data = line_df, ggplot2::aes(x = x, y = y), color = color_map[[prog]], linewidth = cn_line_width)
          }else{
            plot <- plot + 
              ggplot2::geom_path(data = line_df, ggplot2::aes(x = x, y = y), color = color_map[[prog]], linewidth = cn_line_width)
          }
        }else{
          plot <- plot + 
            ggplot2::geom_point(data = point_df, ggplot2::aes(x = x, y = y), size = bin_point_size, color = color_map[[prog]], alpha = bin_point_alpha)
        }
      }
    }
  }
  
  print(plot)
  
  if(output_dir!=FALSE){
    file_path <- paste0(output_dir,"/heal_multisample_genomewide_plot.", device)
    if (file.exists(file_path)) {
      stop(paste0("Already a file at ", file_path,". Exiting.."))
    } else {
      ggplot2::ggsave(filename = file_path, plot = plot, device = device, ...)
      cat(paste0("Plot saved at ", file_path,"."))
    }
  }
  
  if(return_plot==TRUE){
    return(plot)
  }
}

#' Plot a heat map of all CN 
#'
#' @param heal_list List in heal format with CN information (such as output from get_copy_number()).
#' @param prog_ploidy Ploidy of the progenitors (Assumed to be equal. '2' by default).
#' @param sample_label_size Size of the sample labels ('10' by default).
#' @param chr_limits Show chromosome edges with dashed lines ('FALSE' by default).
#' @param chr_limit_thickness Thickness of chromosome limit lines ('0.5' by default).
#' @param chr_labels Show chromosome labels ('FALSE' by default).
#' @param chr_label_size Size of chromosome labels ('4' by default).
#' @param subgenome_limits Show subgenome limits with solid lines ('FALSE' by default).
#' @param subgenome_limit_thickness Thickness of subgenome limit lines ('0.8' by default).
#' @param subgenome_labels Show subgenome labels ('FALSE' by default).
#' @param subgenome_label_size Size of subgenome labels ('5' by default).
#' @param separate_subgenome_plots Plot each subgenome separately ('FALSE' by default).
#'
#' @returns
#' @export
#'
#' @examples
plot_cn_heat <- function(heal_list, prog_ploidy = 2, sample_label_size = 10,
                         chr_limits = FALSE, chr_limit_thickness = 0.5,
                         chr_labels = FALSE, chr_label_size = 4,
                         subgenome_limits = FALSE, subgenome_limit_thickness = 0.8,
                         subgenome_labels = FALSE, subgenome_label_size = 5, 
                         separate_subgenome_plots = FALSE){

  # Check if CN exists
  cn_is_null <- is.null(unlist(lapply(heal_list, function(list) {
    list$CN
  })))
  if (cn_is_null) {
    stop("No CN data table found in heal_list. Make sure the input heal_list contains CN data. \n")
  }

  progenitors <- names(heal_list)
  
  smp_stats <- get_sample_stats(heal_list, sample_type = TRUE)
  polyploid_samples <- smp_stats$sample[smp_stats$type=="polyploid"]
  
  # get CN, chr and prog for each samples. 
  dt_prog_list <- foreach::foreach(prog = progenitors)%do%{
    
    col_keep <- c("chr", polyploid_samples)
    dt_prog <- heal_list[[prog]]$CN[, ..col_keep]
    dt_prog$prog <- rep(prog, nrow(dt_prog))
    return(dt_prog)
  }
  dt_all <- rbindlist(dt_prog_list)
  
  dt_all$bin_index <- 1:nrow(dt_all)
  
  # Make in right format for heat map
  dt_long <- as.data.table(tidyr::pivot_longer(
    data = dt_all,
    cols = -c(bin_index, chr, prog),
    names_to = "sample",
    values_to = "CN"
  ))
  
  # Get position of chr labels 
  rle_chr <- rle(dt_all$chr) 
  end_pos <- c(0, cumsum(rle_chr$lengths))
  chr_lab_pos <- end_pos + (rle_chr$lengths/2)
  
  ### Set size at bottom based on chr label length
  k <- 0.035 # k is a small constant converting mm to data units.
  max_chr_chars <- max(nchar(rle_chr$values))
  extra_bottom_space <- max_chr_chars * chr_label_size * k  
  
  chr_labels <- data.frame(label = rle_chr$values,
                           x_pos = chr_lab_pos[1:length(chr_lab_pos)-1],
                           y_pos = rep( (- extra_bottom_space - 0.5)  / 2, length(rle_chr$values)))
  
  # Get position of chr edges dashed lines
  grad_pos_dt <- data.table::data.table(x = end_pos, xend = end_pos, y = rep(- extra_bottom_space - 1, length(end_pos)), yend = rep(length(polyploid_samples) + 0.5, length(end_pos)))
  
  # Get position of subgenome edges
  sub_g_edge_dt <- data.table::rbindlist(foreach::foreach(prog = progenitors)%do%{
    
    x_pos <- max(dt_all$bin_index[dt_all$prog == prog])
    y_pos <- length(polyploid_samples) + 1.5
    return(data.table::data.table(x = x_pos, xend = x_pos, y = - extra_bottom_space - 1, yend = y_pos))
    
  })
  sub_g_edge_dt <- rbind(data.table::data.table(x = 0, xend = 0, y = - extra_bottom_space - 1, yend = length(polyploid_samples) + 1.5), sub_g_edge_dt)
  
  # Get position of progenitor/subgenome label
  pos_and_prog_dt <- data.table::rbindlist(foreach::foreach(prog = progenitors)%do%{
    
    x_pos <- mean(dt_all$bin_index[dt_all$prog==prog])
    y_pos <- length(polyploid_samples) + 1
    return(data.table::data.table(label=prog, x_pos=x_pos, y_pos=y_pos))
    
  })
  
  # Cap the subgenome labels
  cap_of_plot <- data.table::data.table(x = c(0,0),
                                        xend = c(max(dt_all$bin_index), max(dt_all$bin_index)),
                                        y = c(length(polyploid_samples) + 1.5, - extra_bottom_space - 1),
                                        yend = c(length(polyploid_samples) + 1.5, - extra_bottom_space - 1))
  
  # Make the right plot
  outplot <- ggplot2::ggplot(dt_long, ggplot2::aes(x = bin_index, y = sample, fill = CN)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 2,
      limits = c(0, length(progenitors) * prog_ploidy),
      oob = scales::squish
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = sample_label_size)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::scale_y_discrete(
      expand = ggplot2::expansion(add = c(extra_bottom_space + 1, 0))
    ) +
    ggplot2::geom_text(
      data = chr_labels,
      ggplot2::aes(x = x_pos, y = y_pos, label = label),
      angle = 90,
      size = chr_label_size,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = pos_and_prog_dt,
      ggplot2::aes(x = x_pos, y = y_pos, label = label),
      size = subgenome_label_size,
      inherit.aes = FALSE
    ) +
    ggplot2::labs(
      x = "",
      y = "",
      fill = "Copy Number"
    ) +
    ggplot2::geom_segment(
      data = grad_pos_dt,
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
      linetype = "dashed",
      color = "black",
      size = chr_limit_thickness,
      inherit.aes = FALSE  
    ) +
    ggplot2::geom_segment(
      data = sub_g_edge_dt,
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
      linetype = "solid",
      color = "black",
      size = subgenome_limit_thickness,
      inherit.aes = FALSE  
    ) + 
    ggplot2::geom_segment(
      data = cap_of_plot,
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
      linetype = "solid",
      color = "black",
      size = subgenome_limit_thickness,
      inherit.aes = FALSE  
    )
  
  return(outplot)
}



#' Plotting function for heal_list.
#'
#' @param heal_list List in heal format (such as output from count_heal_data()).
#' @param view_sample The name of a sample to plot (as character)('FALSE' by default).
#' @param output_dir The name of a directory to write all plots to. Will create one if nonexistent.
#' @param n_threads Number of threads to use ('1' by default).
#' @param prog_ploidy Ploidy of the progenitors (Assumed to be equal. '2' by default).
#' @param plot_cn Logical: plot a line indicating infered copy number ('FALSE' by default in CN has been estimated).
#' @param add_bins Logical: plot counts for each bin ('TRUE' by default; normalized in plot_cn=TRUE).
#' @param color_map A vector of colors for each progenitor. If "FALSE" the colors are choosen using viridis().
#' @param specific_chr A vector of characters indicating which chromosomes to plot (plots all by default).
#' @param return_list Logical: return a list of plots if view_sample.
#' @param method Which method was used to assign a copy number to each segment (if plot_cn=TRUE and add_bins!=FALSE). It should be the same method as the one used in get_copy_number() ('global', 'local' or 'manual'. 'global' by default). 
#' @param average The method to compute average to normalize counts (if plot_cn=TRUE and add_bins!=FALSE). It should be the same method as the one used in get_copy_number() ('median' or 'mean'. 'median' by default). 
#' @param average_list Same as in get_copy_number(). 
#' @param ylim_max Maximum y axis value when plotting copy number. 
#' @param device Plot device (as argument to ggplot2::ggsave) ('png' by default).
#' @param width Width of the plot output (as argument to ggplot2::ggsave) ('6' by default).
#' @param height Height of the plot output (as argument to ggplot2::ggsave) ('4' by default).
#' @param units Units of the plot output (as argument to ggplot2::ggsave) ('in' by default).
#' @param ... Any arguments you wish to pass to ggplot2::geom_point()
#'
#' @return Either nothing or a list of plots.
#' @export
#'
plot_bins <- function(heal_list, view_sample = FALSE, output_dir = FALSE,
                      n_threads = 1, prog_ploidy = 2, plot_cn = FALSE, 
                      add_bins = TRUE, add_DNAcopy = FALSE, color_map = FALSE, 
                      specific_chr = FALSE, return_list = FALSE, method = "global", 
                      average = "median", average_list = FALSE, linewidth=2, ylim_max=8, 
                      device = "png", width = 6, height = 4, units = "in", ...) {
  
  # Check input to see if it's appropriate
  if(!is.logical(add_bins)){
    stop("'add_bins' must be either TRUE or FALSE (logical). Exiting..")
  }
  
  if(!is.logical(add_DNAcopy)){
    stop("'add_DNAcopy' must be either TRUE or FALSE (logical). Exiting..")
  }
  
  cn_exist <- unlist(lapply(heal_list, function(list) {
    list$CN
  }))
  if (is.null(cn_exist) && plot_cn == TRUE) {
    stop("ERROR: 'plot_cn' set to 'TRUE' but no CN data table found in heal_list. Setting 'plot_cn' to 'FALSE'. \n")
    plot_cn <- FALSE
  }

  if(add_DNAcopy == TRUE){
    DNAcopy_exist <- unlist(lapply(heal_list, function(list) {
      list$DNAcopy
    }))
    if (is.null(DNAcopy_exist) && add_DNAcopy == TRUE) {
      stop("ERROR: 'add_DNAcopy' set to 'TRUE' but no DNAcopy output found in heal_list. Setting 'add_DNAcopy' to 'FALSE'. \n")
      add_DNAcopy <- FALSE
    }
  }
  
  # Get the averages for normalization
  if(method == "global"){
    average_list <- get_sample_stats(heal_list, method = average)
    sample_names <- names(average_list)
    
  }else if(method == "local"){
    average_list <- get_sample_stats(heal_list, method = paste0("local_", average))
    sample_names <- names(average_list)
    
  }else if (method == "manual" & average_list != FALSE){
    sample_names <- names(average_list)
  }
  
  progenitors <- names(heal_list)

  # Check for manual colour input
  if (isFALSE(color_map)) {
    color_map <- viridis::viridis(length(progenitors))
  } else if (length(color_map) != length(progenitors)) {
    stop("Custom color_map is not of correct length. It should match the number of subgenomes. Exiting..")
  }
  names(color_map) <- progenitors

  # sample type
  smp_type_map <- get_sample_stats(heal_list, sample_type = TRUE)

  # Check if ploting only one sample or plotting all to directory.
  if (view_sample != FALSE) {
    if (sum(sample_names == view_sample) == 0) {
      stop("Sample name not recognized for viewing of counts or coverage. Exiting..")
    }

    if (output_dir == FALSE) {
      cat(paste0("Plotting for ", view_sample, ". \n"))
    } else {
      cat(paste0("Saving ", view_sample, "to ", output_dir, ".", "\n"))
    }
    samples <- view_sample

    smp_type_map <- smp_type_map[smp_type_map$sample == view_sample, ]
  } else if (output_dir == FALSE) {
    stop("No output directory and no 'view_sample' set. One must be set.")
  } else {
    cat(paste0("Saving all samples and chromosomes to ", output_dir, ".", "\n"))
    samples <- sample_names
  }

  # Plot each sample
  for (smp in samples) {
    
    if (output_dir != FALSE) {
      ref_dir <- paste0(output_dir, "/", smp, "/")
      dir.create(ref_dir, showWarnings = FALSE, recursive = TRUE)
    }

    if (smp_type_map$type[smp_type_map$sample == smp] != "polyploid") {
      current_prog <- smp_type_map$type[smp_type_map$sample == smp]
    } else {
      current_prog <- progenitors
    }

    plot_prog_list <- foreach::foreach(prog = current_prog) %do% {
      
      if (output_dir != FALSE) {
        dir.create(paste0(ref_dir, "/", prog, "/"), showWarnings = FALSE, recursive = TRUE)
      }

      if (sum(specific_chr != FALSE) != FALSE) {
        chromo <- intersect(specific_chr, unique(heal_list[[prog]]$bins$chr))
      } else {
        chromo <- unique(heal_list[[prog]]$bins$chr)
      }

      doParallel::registerDoParallel(n_threads)
      plot_list <- foreach::foreach(chr = chromo) %dopar% {
        
        if (output_dir != FALSE) {
          out_file <- paste0(ref_dir, "/", prog, "/", chr, ".", device)
        }

        
        which_bin_rows <- heal_list[[prog]]$bins$chr == chr
        which_cn_rows <- heal_list[[prog]]$CN$chr == chr

        x <- heal_list[[prog]]$bins$start[which_bin_rows]
        y_vec_pts <- heal_list[[prog]]$bins[[smp]][which_bin_rows]

        if(add_DNAcopy == TRUE){
          
          which_index <- heal_list[[prog]]$DNAcopy[[smp]]$output$chrom==chr
          start_vec <- heal_list[[prog]]$DNAcopy[[smp]]$output$loc.start[which_index]
          end_vec <- heal_list[[prog]]$DNAcopy[[smp]]$output$loc.end[which_index]
          
          if(plot_cn == TRUE){
            normalized_mean_vec <- heal_list[[prog]]$DNAcopy[[smp]]$output$seg.mean[which_index] / average_list[[smp]][[prog]] * prog_ploidy
            DNAcopy_mean <- rep(normalized_mean_vec)
          }else{
            DNAcopy_mean <- rep(heal_list[[prog]]$DNAcopy[[smp]]$output$seg.mean[which_index])
          }
          
          DNAcopy_df <- data.frame(x_start = start_vec, x_end = end_vec, mean = DNAcopy_mean, color_legend = rep("Segment Mean", length(end_vec)))
          color_map <- c(color_map, "Segment Mean" = "red")
        }
        
        if (plot_cn == TRUE) {
          
          y_vec_line <- heal_list[[prog]]$CN[[smp]][which_cn_rows]
          
          if (add_bins == TRUE) {
          
            y_vec_pts <- (y_vec_pts / average_list[[smp]][[prog]]) * prog_ploidy
            
            color_map <- c(color_map, "Normalized Count" = "black")
             
          }

          plot_df <- data.frame(start = x, counts = y_vec_pts, copy = y_vec_line, Legend = rep("Copy Number", length(y_vec_line)))

          final_color_map <- color_map
          names(final_color_map)[names(color_map)==prog] <- "Copy Number"
          
          bin_plot <- ggplot2::ggplot() +
            ggplot2::geom_line(data = plot_df, ggplot2::aes(x = x, y = copy, color = Legend), linewidth = linewidth) +
            ggplot2::geom_point(data = plot_df, ggplot2::aes(x = x, y = counts, color = Legend), ..., size = 1, alpha = 0.1) +
            ggplot2::theme_minimal() +
            ggplot2::scale_color_manual(values = final_color_map) +
            ggplot2::ylim(0, ylim_max) +
            ggplot2::labs(title = paste0(chr, "; ", smp, "; ", prog), x = "Position", y = "Copy Number") +
            ggplot2::theme_bw()
        
          if(add_bins == TRUE){
            bin_plot <- bin_plot +
              ggplot2::geom_point(ggplot2::aes(x = 0, y = -100, color = "Normalized Count"), shape = 16) 
          }
          
          if(add_DNAcopy == TRUE){
            bin_plot <- bin_plot + 
              ggplot2::geom_segment(
                data = DNAcopy_df,
                ggplot2::aes(x = x_start, xend = x_end, y = mean, yend = mean, color = color_legend),
                linewidth = 0.7) +
              ggplot2::scale_color_manual(values = final_color_map)
          }

          if (output_dir != FALSE) {
            ggplot2::ggsave(filename = out_file, bin_plot, device = device, width = width, height = height, units = units)
          } else {
            return(bin_plot)
          }
          
        } else {
          
          plot_df <- data.frame(start = x, counts = y_vec_pts, progenitor = rep(prog, length(y_vec_pts)))

          bin_plot <- ggplot2::ggplot() +
            ggplot2::geom_point(data = plot_df, ggplot2::aes(x = x, y = counts, color = progenitor), ..., size = 1, alpha = 1) +
            ggplot2::theme_minimal() +
            ggplot2::scale_color_manual(values = color_map) +
            ggplot2::labs(title = paste(chr, smp), x = "Position", y = "Counts") +
            ggplot2::theme_bw()
          
          if(add_DNAcopy == TRUE){
            bin_plot <- bin_plot + 
              ggplot2::geom_segment(
                data = DNAcopy_df,
                ggplot2::aes(x = x_start, xend = x_end, y = mean, yend = mean, color = color_legend),
                linewidth = 0.7) +
              ggplot2::scale_color_manual(values = color_map)
          }

          if (output_dir != FALSE) {
            ggplot2::ggsave(filename = out_file, bin_plot, device = device, width = width, height = height, units = units)
          } else {
            return(bin_plot)
          }
        }
      }
      doParallel::stopImplicitCluster()
      
      if (view_sample != FALSE) {
        if (return_list == TRUE) {
          names(plot_list) <- chromo
          return(plot_list)
        } else {
          lapply(plot_list, print)
        }
      }
    }
  }
  if (view_sample != FALSE) {
    if (return_list == TRUE) {
      names(plot_prog_list) <- current_prog
      return(plot_prog_list)
    }
  }
}

#' Plotting function for heal alignment
#'
#' @param heal_list Output list in heal format (such as output from count_heal_data()).
#' @param alignment A heal alignment object created with get_heal_alignment().
#' @param view_sample The name of a sample to plot (as character)('FALSE' by default). If output_dir!=FALSE, the plots will be saved, otherwise they are simply shown. 
#' @param output_dir The name of a directory to write all plots to. Will create one if nonexistent.
#' @param n_threads Number of threads to use ('1' by default).
#' @param add_bins Add points for bins ('FALSE' by default). If "ref" the bins normalized copy number (divided by average (median by default)). If "all" the bins overlapping with anchors are also added at the starting position of the reference anchor. If "FALSE" no bins are added.
#' @param method Which method was used to assign a copy number to each segment (if plot_cn=TRUE and add_bins!=FALSE). It should be the same method as the one used in get_copy_number() ('global', 'local' or 'manual'. 'global' by default). 
#' @param average The method to compute average to normalize counts (if add_bins!=FALSE). It should be the same method as the one used in get_copy_number() ('median' or 'mean'. 'median' by default). 
#' @param average_list Same as in get_copy_number(). 
#' @param prog_ploidy Ploidy of the progenitors (Assumed to be equal. '2' by default).
#' @param color_map A vector of colors for each progenitor. If "FALSE" the colors are choosen using viridis().
#' @param specific_chr A vector of characters indicating which chromosomes to plot (plots all by default).
#' @param return_list Logical: return a list of plots if view_sample.
#' @param ylim_max Maximum y axis value when plotting copy number. 
#' @param device Plot device (as argument to ggplot2::ggsave) ('png' by default).
#' @param width Width of the plot output (as argument to ggplot2::ggsave) ('6' by default).
#' @param height Height of the plot output (as argument to ggplot2::ggsave) ('4' by default).
#' @param units Units of the plot output (as argument to ggplot2::ggsave) ('in' by default).
#' @param ... Any arguments you wish to pass to ggplot2::geom_point()
#'
#' @return Either nothing or a list of plots.
#' @export
#'
plot_alignment <- function(heal_list, alignment, method = "global",
                           average = "median", average_list=FALSE, view_sample = FALSE,
                           output_dir = FALSE, n_threads = 1, add_bins = FALSE,
                           prog_ploidy = 2, color_map = FALSE, specific_chr = FALSE,
                           return_list = FALSE, ylim_max=8, 
                           device = "png", width = 6, height = 4, units = "in", ...) {
  
  if (!add_bins %in% c(FALSE, "ref", "all")) {
    stop("Please input a valid 'add_bins' value. Allowed are: FALSE, 'ref' and 'all'. Exiting..")
  }

  polyploid_samples <- names(alignment)

  # Define averages depending on the method set
  if(method == "global"){
    average_list <- get_sample_stats(heal_list, method = average)
    sample_names <- names(average_list)
    
  }else if(method == "local"){
    average_list <- get_sample_stats(heal_list, method = paste0("local_", average))
    sample_names <- names(average_list)
    
  }else if (method == "manual" & average_list != FALSE){
    sample_names <- names(average_list)
  }
  
  progenitors <- names(heal_list)

  # Define colours: if not set manually, sample with viridis
  if (isFALSE(color_map)) {
    color_map <- viridis::viridis(length(progenitors))
  } else if (length(color_map) != length(progenitors)) {
    stop("Custom color_map is not of correct length. It should match the number of subgenomes. Exiting..")
  }
  names(color_map) <- progenitors
  
  # Parse plotting options and make sure they are compatible with one another.
  if (view_sample != FALSE) {
    if (sum(polyploid_samples == view_sample) == 0) {
      stop("Sample name not recognized (or not polyploid) for viewing of alignment. Exiting..")
    }

    if (output_dir == FALSE) {
      cat(paste0("Plotting for ", view_sample, ". \n"))
    } else {
      cat(paste0("Saving ", view_sample, "to ", output_dir, ".", "\n"))
    }

    polyploid_samples <- view_sample
  } else if (output_dir == FALSE) {
    stop("No output directory and no 'view_sample' set. One must be set. Exiting..")
  } else {
    cat(paste0("Plotting all samples and chromosomes to ", output_dir, ".", "\n"))
  }

  # Define offset for line visibility
  n <- length(progenitors)
  step <- 0.1
  if (n %% 2 == 1) {
    # Odd length: symmetric with 0 in the center
    offsets <- seq(from = -step * floor(n / 2), to = step * floor(n / 2), by = step)
  } else {
    # Even length: symmetric around 0, but no 0 exactly
    half <- n / 2
    offsets <- c(seq(from = -step * half + step / 2, by = step, length.out = half),
      seq(from = step / 2, by = step, length.out = half))
  }
  names(offsets) <- progenitors
  
  # Make plot for each sample
  foreach::foreach(smp = polyploid_samples) %do% {
    
    # We progressively build the outputs and paths if an output dir is set. 
    if (output_dir != FALSE) {
      smp_dir <- paste0(output_dir, "/", smp, "/")
      dir.create(smp_dir, showWarnings = FALSE, recursive = TRUE)
    }

    # Plot along every progenitor gene order successively.
    plot_ref_list <- foreach::foreach(ref = progenitors) %do% {
      
      alt_gnms <- setdiff(progenitors, ref)

      if (output_dir != FALSE) {
        ref_dir <- paste0(smp_dir, "/", ref, "/")
        dir.create(ref_dir, showWarnings = FALSE, recursive = TRUE)
      }

      if (sum(specific_chr != FALSE) != FALSE) {
        chromo <- intersect(specific_chr, unique(heal_list[[ref]]$bins$chr))
      } else {
        chromo <- unique(heal_list[[ref]]$bins$chr)
      }

      doParallel::registerDoParallel(n_threads)
      plot_list <- foreach::foreach(chr = chromo) %dopar% {

        if (output_dir != FALSE) {
          out_file <- paste0(ref_dir, "/", chr, ".", device)
        }

        # Get the current reference line data
        ref_chr_col_name <- paste0("chr_", ref)
        ref_start_col_name <- paste0("start_", ref)
        ref_end_col_name <- paste0("end_", ref)
        ref_cn_col_name <- paste0("cn_", ref)

        which_rows_aln_dt <- alignment[[smp]][[ref_chr_col_name]] == chr

        x_line <- (alignment[[smp]][[ref_start_col_name]][which_rows_aln_dt] + alignment[[smp]][[ref_end_col_name]][which_rows_aln_dt]) / 2
        x_vec_line <- x_line 
        y_vec_line <- alignment[[smp]][[ref_cn_col_name]][which_rows_aln_dt] + offsets[ref]
        subgnm_group <- rep(ref, length(x_vec_line))

        # Add the alt line info
        for (alt in alt_gnms) {
          alt_cn_col_name <- paste0("cn_", alt)
          y_alt <- alignment[[smp]][[alt_cn_col_name]][which_rows_aln_dt] + offsets[alt]

          subgnm_group <- c(subgnm_group, rep(alt, length(y_alt)))
          x_vec_line <- c(x_vec_line, x_line) # same coordinates
          y_vec_line <- c(y_vec_line, y_alt) 
        }

        lines_df <- data.frame(start = x_vec_line, copy = y_vec_line, Legend = subgnm_group)
        
        # if add bins, create a data table with bin information
        if (add_bins != FALSE) {
          which_rows_bins_dt <- heal_list[[ref]]$bins$chr == chr

          x_pts <- (heal_list[[ref]]$bins$start[which_rows_bins_dt] + heal_list[[ref]]$bins$end[which_rows_bins_dt]) / 2
          x_vec_pts <- x_pts
          y_vec_pts <- heal_list[[ref]]$bins[[smp]][which_rows_bins_dt]
          
          # normalize
          y_vec_pts <- (y_vec_pts / average_list[[smp]][[ref]]) * prog_ploidy 
          
          subgnm_group <- rep(ref, length(y_vec_pts))
          
          color_map <- c(color_map, "Normalized Count" = "black")
          
          # add bins of alt if all
          if (add_bins == "all") {
            for (alt in alt_gnms) {
              alt_bin_col_name <- paste0("bin_index_", alt)
              
              # Only if the current ref chromosome has synteny to alt subgenome
              if(nrow(alignment[[smp]][which_rows_aln_dt, ]) != 0){
                
                # Get the alt bin index and the corresponding x position on ref. 
                cn_bindex_dt_list <- apply(alignment[[smp]][which_rows_aln_dt, ], 1, function(row) {
                  bindex_vec <- as.numeric(unlist(strsplit(row[[alt_bin_col_name]], ",")))
                  midpoint <- (as.numeric(row[[ref_start_col_name]]) + as.numeric(row[[ref_end_col_name]])) / 2
                  start_vec <- rep(midpoint, length(bindex_vec))
                  return(data.table::data.table(bin_index = bindex_vec, ref_start = start_vec))
                })
  
                bin_cn_dt <- unique(data.table::rbindlist(cn_bindex_dt_list))
                
                # Here we get the bin count of the desired alt bins by
                # first getting the chr and start position in the CN data table. 
                merge_dt <- data.table::data.table(chr = heal_list[[alt]]$CN$chr[bin_cn_dt$bin_index], start = heal_list[[alt]]$CN$start[bin_cn_dt$bin_index])
                # then merging this to the bins data table. 
                merge_dt <- merge(merge_dt, heal_list[[alt]]$bins, sort = FALSE)
  
                x_pts_alt <- as.numeric(bin_cn_dt$ref_start)
                y_pts_alt <- merge_dt[[smp]]
                
                y_pts_alt <- (y_pts_alt / average_list[[smp]][[alt]]) * prog_ploidy
                
                x_vec_pts <- c(x_vec_pts, x_pts_alt)
                y_vec_pts <- c(y_vec_pts, y_pts_alt)
                subgnm_group <- c(subgnm_group, rep(alt, length(y_pts_alt)))
              }
            }
          }

          pts_df <- data.frame(start = x_vec_pts, points = y_vec_pts, Legend = subgnm_group)

          bin_plot <- ggplot2::ggplot() +
            ggplot2::geom_line(data = lines_df, ggplot2::aes(x = start, y = copy, color = Legend), linewidth = 2) +
            ggplot2::geom_point(data = pts_df, ggplot2::aes(x = start, y = points, color = Legend), ..., size = 1, alpha = 0.1) +
            ggplot2::geom_point(ggplot2::aes(x = 0, y = -100, color = "Normalized Count"), shape = 16) +
            ggplot2::theme_minimal() +
            ggplot2::scale_color_manual(values = color_map) +
            ggplot2::ylim(-0.1, ylim_max) +
            ggplot2::labs(title = paste(chr, "; ", smp, "; ", ref), x = "Position", y = "Copy Number") +
            ggplot2::theme_bw()

          if (output_dir != FALSE) {
            ggplot2::ggsave(filename = out_file, bin_plot, device = device, width = width, height = height, units = units)
          } else {
            return(bin_plot)
          }
        } else {
          aln_plot <- ggplot2::ggplot() +
            ggplot2::geom_line(data = lines_df, ggplot2::aes(x = start, y = copy, color = Legend), linewidth = 2) +
            ggplot2::theme_minimal() +
            ggplot2::scale_color_manual(values = color_map) +
            ggplot2::ylim(-0.1, ylim_max) +
            ggplot2::labs(title = paste(chr, "; ", smp, "; ", ref), x = "Position", y = "Copy Number") +
            ggplot2::theme_bw()

          if (output_dir != FALSE) {
            ggplot2::ggsave(filename = out_file, aln_plot, device = device, width = width, height = height, units = units)
          } else {
            return(aln_plot)
          }
        }
      }
      doParallel::stopImplicitCluster()

      if (view_sample != FALSE) {
        if (return_list == TRUE) {
          names(plot_list) <- chromo
          return(plot_list)
        } else {
          lapply(plot_list, print)
        }
      }
    }
  }
  if (view_sample != FALSE) {
    if (return_list == TRUE) {
      names(plot_ref_list) <- progenitors
      return(plot_ref_list)
    }
  }
}


#' Plot empirical concordant densities
#'
#' @param densities A list of densities for each polyploid sample (output from get_concordant_density()).
#' @param view_sample The name of a sample to plot (as character)('FALSE' by default).
#' @param output_dir The name of a directory to write all plots to. Will create one if nonexistent.
#' @param show_discordant Logical defining if bins involved in a discordant anchor set should be added on the plot.
#' @param heal_list Output list in heal format (such as output from count_heal_data()).
#' @param alignment A heal alignment object created with get_heal_alignment().
#' @param corrected_alignment A density corrected alignment ('FALSE' by default). If one is given and show_discordant==TRUE then the new state of corrected points is shown above the previous state.
#' @param ylim_max Maximum ylim value ('FALSE' by default). If FALSE then the value is set for each sample as 3 times the standard deviation above the mean.
#' @param color_vec A vector of colors for each copy number class from 0 to the ploidy of the progenitors times the number of progenitors ('FALSE' by default). If "FALSE" the colors are choosen using viridis().
#' @param prog_ploidy Ploidy of the progenitors (Assumed to be equal. '2' by default)
#' @param n_threads Number of threads to use ('1' by default).
#' @param normalize Use normalize count values (either FALSE, 'local' or 'manual'. FALSE by default). The 'local' and 'manual' entries have the same meaning as for 'method' in 'get_copy_number()'.  
#' @param average Average measure used normalize each bin count ('median' or 'mean'. 'median' by default).
#' @param average_list A list with one element per sample with each containing one element per subgenome with values used for normalization (same as in 'get_copy_number()').

#'
#' @return Nothing. Plots are shown and/or saved to output_dir.
#' @export
#'
#' @importFrom data.table :=
plot_densities <- function(densities, view_sample = FALSE, output_dir = FALSE,
                           show_discordant = FALSE, heal_list = FALSE, alignment = FALSE,
                           corrected_alignment = FALSE, ylim_max = FALSE, color_vec = FALSE,
                           prog_ploidy = 2, n_threads = 1, normalize = FALSE,
                           average = "median", average_list = FALSE){
  
  is_align_and_count_data <- is.list(alignment) & is.list(heal_list)
  if (show_discordant == TRUE & !is_align_and_count_data) {
    stop("show_discordant is TRUE but no valid heal_list or aligment provided. Exiting..")
  }

  polyploid_samples <- names(densities)
  progenitors <- names(heal_list)

  
  if(normalize == "local"){
    local_averages <- get_sample_stats(heal_list, method = paste0("local_", average))
    
  }else if (normalize == "manual" & average_list != FALSE){
    if (!identical(names(average_list), polyploid_samples)){
      stop("Invalid average_list input. When method == 'manual' you must provide an average_list with values of median for each sample subgenome.")
    }
  }else if (normalize == "manual" & average_list == FALSE){
    stop("No 'average_list' provided. When method == 'manual' you must provide an 'average_list' with values of median for each sample subgenome.")
  }
  
  
  if (view_sample != FALSE) {
    n_threads <- 1
    if (sum(polyploid_samples == view_sample) == 0) {
      stop("Sample name not recognized (or not polyploid) for viewing of alignment. Exiting..")
    }

    cat(paste0("Plotting for ", view_sample, ". \n"))
    polyploid_samples <- view_sample
  } else if (output_dir == FALSE) {
    stop("No output directory and no 'view_sample' set. One must be set. Exiting..")
  } else {
    cat(paste0("Plotting all samples and chromosomes to ", output_dir, ".", "\n"))
  }

  doParallel::registerDoParallel(n_threads)
  catch_all <- foreach::foreach(smp = polyploid_samples) %dopar% {
    if (output_dir != FALSE) {
      dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
      output_file <- paste0(output_dir, "/", smp, "_density.png")
      grDevices::png(output_file)
    }

    if (ylim_max == FALSE) {
      if (is.list(heal_list) == FALSE) {
        
        stop("No ylim_max and no heal_list. Set at least one. Exiting..")
        
      } else {
        
        counts_list <- foreach::foreach(prog = progenitors)%do%{
          if(normalize == FALSE){
            counts <- heal_list[[prog]]$bins[[smp]] 
          } else if(normalize == "local"){
            counts <- heal_list[[prog]]$bins[[smp]]  / local_averages[[smp]][[prog]] * prog_ploidy
          }else if(normalize == "manual"){
            counts <- heal_list[[prog]]$bins[[smp]]  / average_list[[smp]][[prog]] * prog_ploidy
          }
          return(counts)
        }
        counts_vec <- unlist(counts_list)
        
        ylim_max <- mean(counts_vec, na.rm = TRUE) + 3 * stats::sd(counts_vec, na.rm = TRUE)
      }
    }

    cn_labels <- names(densities[[smp]])

    if (isFALSE(color_vec)) {
      color_vec <- viridis::viridis(n = length(cn_labels))
    } else {
      if (length(cn_labels) != length(color_vec)) {
        stop("Color vector length not matching number of copy number categories. Exiting..")
      }
    }

    names(color_vec) <- cn_labels
    
    if(normalize == FALSE){
      graphics::contour(densities[[smp]][[cn_labels[1]]], ylim = c(-5, ylim_max), col = color_vec[cn_labels[1]], main = smp)
    }else{
      graphics::contour(densities[[smp]][[cn_labels[1]]], ylim = c(0, ylim_max), col = color_vec[cn_labels[1]], main = smp)
    }
    
    if(length(cn_labels)>1){
      for (i in 2:length(cn_labels)) {
        graphics::contour(densities[[smp]][[cn_labels[i]]], ylim = c(0, ylim_max), col = color_vec[cn_labels[i]], add = TRUE)
      }
    }

    if (show_discordant == TRUE) {
      which_discord <- alignment[[smp]]$status == "discordant"

      per_prog_dt_list <- foreach::foreach(prog = progenitors) %do% {
        cn_col_name <- paste0("cn_", prog)
        bin_col_name <- paste0("bin_index_", prog)
        col_keep <- c(cn_col_name, bin_col_name)

        cn_bindex_dt_list <- apply(alignment[[smp]][which_discord, ..col_keep], 1, function(row) {
          bindex_vec <- as.numeric(unlist(strsplit(row[[bin_col_name]], ",")))
          cn_vec <- rep(row[[cn_col_name]], length(bindex_vec))
          return(data.table::data.table(bin_index = bindex_vec, cn = cn_vec))
        })

        bin_cn_dt <- unique(data.table::rbindlist(cn_bindex_dt_list))

        merge_dt <- data.table::data.table(chr = heal_list[[prog]]$CN$chr[bin_cn_dt$bin_index], start = heal_list[[prog]]$CN$start[bin_cn_dt$bin_index])

        merge_dt <- merge(heal_list[[prog]]$bins, merge_dt, by = c("chr", "start"))

        col_keep <- c("gc_content", smp)
        
        if(normalize == "local"){
          merge_dt[[smp]] <- merge_dt[[smp]] / local_averages[[smp]][[prog]] * prog_ploidy
          
        }else if(normalize == "manual"){
          merge_dt[[smp]] <- merge_dt[[smp]] / average_list[[smp]][[prog]] * prog_ploidy
          
        }

        return(cbind(merge_dt[, ..col_keep], cn = bin_cn_dt$cn))
      }

      points_dt <- data.table::rbindlist(per_prog_dt_list)

      points_dt$cn <- paste0("cn_", points_dt$cn)

      points_dt[, is_allowed := cn %in% cn_labels]

      color_vec[["other"]] <- "black"
      points_dt$cn[points_dt$is_allowed == FALSE] <- "other"

      graphics::points(points_dt$gc_content, points_dt[[smp]], col = color_vec[points_dt$cn], pch = 16, cex = 0.9)

      if (is.list(corrected_alignment) == TRUE) {
        which_corrected <- grep("corrected_", corrected_alignment[[smp]]$status)

        per_prog_dt_list <- foreach::foreach(prog = progenitors) %do% {
          cn_col_name <- paste0("cn_", prog)
          bin_col_name <- paste0("bin_index_", prog)
          col_keep <- c(cn_col_name, bin_col_name)

          cn_bindex_dt_list <- apply(corrected_alignment[[smp]][which_corrected, ..col_keep], 1, function(row) {
            bindex_vec <- as.numeric(unlist(strsplit(row[[bin_col_name]], ",")))
            cn_vec <- rep(row[[cn_col_name]], length(bindex_vec))
            return(data.table::data.table(bin_index = bindex_vec, cn = cn_vec))
          })

          bin_cn_dt <- unique(data.table::rbindlist(cn_bindex_dt_list))

          merge_dt <- data.table::data.table(chr = heal_list[[prog]]$CN$chr[bin_cn_dt$bin_index], start = heal_list[[prog]]$CN$start[bin_cn_dt$bin_index])
          merge_dt <- merge(heal_list[[prog]]$bins, merge_dt, by = c("chr", "start"))
          merge_dt <- merge(heal_list[[prog]]$CN, merge_dt, by = c("chr", "start", "gc_content"), suffixes = c("_CN", "_count"))
          col_keep <- c("gc_content", paste0(smp, c("_CN", "_count")))
          merge_dt <- merge_dt[, ..col_keep]
          data.table::setnames(merge_dt, old = col_keep, c("gc_content", "cn_original", "count"))
          
          if(normalize == "local"){
            merge_dt$count <- merge_dt$count / local_averages[[smp]][[prog]] * prog_ploidy
            
          }else if(normalize == "manual"){
            merge_dt$count <- merge_dt$count / average_list[[smp]][[prog]] * prog_ploidy
            
          }
          
          
          return(cbind(merge_dt, cn_corrected = bin_cn_dt$cn))
        }

        points_dt <- data.table::rbindlist(per_prog_dt_list)

        points_dt$cn_original <- paste0("cn_", points_dt$cn_original)
        points_dt$cn_corrected <- paste0("cn_", points_dt$cn_corrected)

        points_dt <- points_dt[points_dt$cn_corrected != points_dt$cn_original, ]

        correction_improved_bin_density <- c()
        for (i in 1:nrow(points_dt)) {
          row <- points_dt[i, ]

          original_cn <- row[["cn_original"]]
          if (original_cn %in% cn_labels) {
            x_idx_orig <- which.min(abs(densities[[smp]][[original_cn]]$x - row[["gc_content"]]))
            y_idx_orig <- which.min(abs(densities[[smp]][[original_cn]]$y - row[["count"]]))
            original_density <- densities[[smp]][[original_cn]]$z[x_idx_orig, y_idx_orig]

            corrected_cn <- row[["cn_corrected"]]
            x_idx_cor <- which.min(abs(densities[[smp]][[corrected_cn]]$x - row[["gc_content"]]))
            y_idx_cor <- which.min(abs(densities[[smp]][[corrected_cn]]$y - row[["count"]]))
            corrected_density <- densities[[smp]][[corrected_cn]]$z[x_idx_cor, y_idx_cor]

            correction_improved_bin_density <- c(correction_improved_bin_density, corrected_density > original_density)
          } else {
            correction_improved_bin_density <- c(correction_improved_bin_density, TRUE)
          }
        }

        graphics::points(points_dt$gc_content[correction_improved_bin_density], points_dt$count[correction_improved_bin_density], col = color_vec[points_dt$cn_corrected[correction_improved_bin_density]], pch = 18, cex = 0.6)
      }
    }
    # reinitialize color_vec
    color_vec <- color_vec[names(color_vec) != "other"]

    if (output_dir != FALSE) {
      grDevices::dev.off()
    }
  }
  doParallel::stopImplicitCluster()
}


#' Plot pairwise heatmap
#'
#' @param alignment A heal alignment object created with get_heal_alignment().
#' @param view_samples A vector of sample names to plot (as character)('FALSE' by default).
#' @param output_dir The name of a directory to write all plots to. Will create one if nonexistent.
#' @param xrange A sequence of value to show for the x axis, such as start:finish (example 1:5) ('FALSE' by default). Make sure to set both xrange and yrange.
#' @param yrange A sequence of value to show for the y axis, such as start:finish (example 1:5) ('FALSE' by default). Make sure to set both xrange and yrange.
#' @param device Plot device (as argument to ggplot2::ggsave) ('pdf' by default).
#' @param ... Any arguments you wish to pass to ggplot2::ggsave().
#'
#' @return Nothing. Plots are shown and/or saved to output_dir.
#' @export
#'
plot_heal_heat_map <- function(alignment, view_samples = FALSE,
                               output_dir = FALSE, xrange = FALSE,
                               yrange = FALSE, device = "pdf", ...) {
  
  polyploid_samples <- names(alignment)


  if (!isFALSE(view_samples)) {
    if (length(intersect(polyploid_samples, view_samples)) == 0) {
      stop("Sample names not recognized (or not polyploid) for viewing of alignment. Exiting..")
    }

    cat(paste("Plotting for:", view_samples, " \n"))
    polyploid_samples <- view_samples
  } else if (isFALSE(output_dir)) {
    stop("No output directory and no 'view_samples' set. One must be set. Exiting..")
  } else {
    cat(paste0("Plotting all samples to ", output_dir, ".", "\n"))
  }


  for (smp in polyploid_samples) {
    cn_cols <- grep("cn_", colnames(alignment[[smp]]), value = TRUE)
    all_pairs <- utils::combn(cn_cols, 2)
    
    # create output directory for this samples
    if (output_dir != FALSE) {
      smp_out_dir <- paste0(output_dir, "/", smp)
      dir.create(smp_out_dir, showWarnings = FALSE, recursive = TRUE)
    }
    
    # for all pairs of progenitors
    for (i in 1:ncol(all_pairs)) {
      
      pair <- all_pairs[, i]
      progenitors <- sub("cn_", "", pair)

      input_dt <- data.table::data.table(x = alignment[[smp]][[pair[1]]], y = alignment[[smp]][[pair[2]]])

      count_table <- as.data.frame(table(input_dt))
      colnames(count_table) <- c("x", "y", "count")
     
      # Use the manually provided ranges if they exist. 
      if(!is.logical(xrange)){
        if(!is.logical(yrange)){
          
          xrange <- as.character(xrange)
          yrange <- as.character(yrange)
          
          full_grid <- expand.grid(x = xrange, y = yrange, stringsAsFactors = FALSE)
          
          count_table_full <- merge(full_grid, count_table, by = c("x", "y"), all.x = TRUE)
          
          count_table_full$count[is.na(count_table_full$count)] <- 0
          count_table <- count_table_full
        
        }else{
          
          xrange <- as.character(xrange)
          y_min <- min(as.numeric(count_table$y))
          y_max <- max(as.numeric(count_table$y))
          yrange <- as.character(y_min:y_max)
          
          full_grid <- expand.grid(x = xrange, y = yrange, stringsAsFactors = FALSE)
          
          count_table_full <- merge(full_grid, count_table, by = c("x", "y"), all.x = TRUE)
          
          count_table_full$count[is.na(count_table_full$count)] <- 0
          count_table <- count_table_full 
        }
      }else if(!is.logical(yrange)){
        
        x_min <- min(as.numeric(count_table$x))
        x_max <- max(as.numeric(count_table$x))
        xrange <- as.character(x_min:x_max)
        yrange <- as.character(yrange)
        
        full_grid <- expand.grid(x = xrange, y = yrange, stringsAsFactors = FALSE)
        
        count_table_full <- merge(full_grid, count_table, by = c("x", "y"), all.x = TRUE)
        
        count_table_full$count[is.na(count_table_full$count)] <- 0
        count_table <- count_table_full
      }
      
    
      plot <- ggplot2::ggplot(count_table, ggplot2::aes(x = x, y = y, fill = count)) +
        ggplot2::geom_tile(color = "black") +
        ggplot2::scale_fill_gradient(low = "blue", high = "red", na.value = "grey50") +
        ggplot2::labs(
          title = paste0("Heatmap of Counts for Homoeolog Pairs in ", smp),
          x = paste("Infered CN in", progenitors[1]), y = paste("Infered CN in", progenitors[2])
        ) +
        ggplot2::theme_minimal()
      
      if (output_dir != FALSE) {
        out_file <- paste0(output_dir, "/", smp, "/", pair[1], "_vs_", pair[2], "_heat.", device)
        ggplot2::ggsave(filename = out_file, plot, device = device, width = 6, height = 4, units = "in", ...)
      } else {
        print(plot)
      }
    }
  }
}

utils::globalVariables(c("a", "..col_keep", "cn", "i", "status", "..bin_index_col", "..which_relevant", "..cn_col", "bin_start", "anchor_start", "anchor_end", "gene_id", "ref", "copy", "subgenome", "points", "progenitor", "counts", "count", "is_allowed", "x", "y"))
