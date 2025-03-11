#' Plotting function for heal_list.
#'
#' @param heal_list Output list in heal format (such as output from count_heal_data()).
#' @param view_sample The name of a sample to plot (as character)('FALSE' by default).
#' @param output_dir The name of a directory to write all plots to. Will create one if nonexistent.
#' @param n_threads Number of threads to use ('1' by default).
#' @param prog_ploidy Ploidy of the progenitors (Assumed to be equal. '2' by default).
#' @param plot_cn Logical: plot a line indicating infered copy number ('FALSE' by default in CN has been estimated).
#' @param add_bins Logical: plot counts for each bin ('TRUE' by default; normalized in plot_cn=TRUE).
#' @param color_map A vector of colors for each progenitor. If "FALSE" the colors are choosen using rainbow().
#' @param specific_chr A vector of characters indicating which chromosomes to plot (plots all by default).
#' @param return_list Logical: return a list of plots if view_sample.
#'
#' @return Either nothing or a list of plots.
#' @export
#'
plot_bins <- function(heal_list, view_sample = FALSE, output_dir = FALSE, n_threads = 1, prog_ploidy = 2, plot_cn = FALSE, add_bins = TRUE, color_map = FALSE, specific_chr = FALSE, return_list = FALSE) {
  cn_exist <- unlist(lapply(heal_list, function(list) {
    list$CN
  }))
  if (is.null(cn_exist) && plot_cn == TRUE) {
    cat("ERROR: plot_cn set to 'TRUE' but no CN data table found in heal_list. Setting plot_cn to 'FALSE'.")
    plot_cn <- FALSE
  }

  sample_averages <- get_sample_stats(heal_list)
  progenitors <- names(heal_list)

  if (isFALSE(color_map)) {
    color_map <- grDevices::rainbow(length(progenitors), s = 0.7)
  } else if (length(color_map) != length(progenitors)) {
    cat("ERROR: Custom color_map is not of correct length. It should match the number of subgenomes. Exiting..")
    return()
  }
  names(color_map) <- progenitors

  smp_type_map <- get_sample_stats(heal_list, sample_type = TRUE)

  if (view_sample != FALSE) {
    if (sum(names(sample_averages) == view_sample) == 0) {
      cat("ERROR: Sample name not recognized for viewing of alignment. Exiting..")
      return()
    }

    if (output_dir == FALSE) {
      cat(paste0("Plotting for ", view_sample, ". \n"))
    } else {
      cat(paste0("Saving ", view_sample, "to ", output_dir, "."))
    }
    samples <- view_sample

    smp_type_map <- smp_type_map[smp_type_map$sample == view_sample, ]
  } else if (output_dir == FALSE) {
    cat("ERROR: no output directory and no 'view_sample' set. One must be set.")
    return()
  } else {
    cat(paste0("Saving all samples and chromosomes to ", output_dir, "."))
    samples <- names(sample_averages)
  }

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
          out_file <- paste0(ref_dir, "/", prog, "/", chr, ".png")
        }

        
        which_bin_rows <- heal_list[[prog]]$bins$chr == chr
        which_cn_rows <- heal_list[[prog]]$CN$chr == chr

        x <- heal_list[[prog]]$bins$start[which_bin_rows]
        y_vec_pts <- heal_list[[prog]]$bins[[smp]][which_bin_rows]

        if (plot_cn == TRUE) {
          y_vec_line <- heal_list[[prog]]$CN[[smp]][which_cn_rows]
          if (add_bins == TRUE) {
            y_vec_pts <- (y_vec_pts / sample_averages[[smp]]) * prog_ploidy
          }

          plot_df <- data.frame(start = x, counts = y_vec_pts, copy = y_vec_line, progenitor = rep(prog, length(y_vec_line)))

          bin_plot <- ggplot2::ggplot() +
            ggplot2::geom_line(data = plot_df, ggplot2::aes(x = x, y = copy, color = progenitor), linewidth = 2) +
            ggplot2::geom_point(data = plot_df, ggplot2::aes(x = x, y = counts, color = progenitor), size = 1, alpha = 0.1) +
            ggplot2::theme_minimal() +
            ggplot2::scale_color_manual(values = color_map) +
            ggplot2::ylim(0, 8) +
            ggplot2::labs(title = paste(chr, smp), x = "Position", y = "Copy Number") +
            ggplot2::theme_bw()

          if (output_dir != FALSE) {
            ggplot2::ggsave(filename = out_file, bin_plot, device = "png", width = 6, height = 4, units = "in")
          } else {
            return(bin_plot)
          }
        } else {
          plot_df <- data.frame(start = x, counts = y_vec_pts, progenitor = rep(prog, length(y_vec_pts)))

          bin_plot <- ggplot2::ggplot() +
            ggplot2::geom_point(data = plot_df, ggplot2::aes(x = x, y = counts, color = progenitor), size = 1, alpha = 1) +
            ggplot2::theme_minimal() +
            ggplot2::scale_color_manual(values = color_map) +
            ggplot2::labs(title = paste(chr, smp), x = "Position", y = "Counts") +
            ggplot2::theme_bw()

          if (output_dir != FALSE) {
            ggplot2::ggsave(filename = out_file, bin_plot, device = "png", width = 6, height = 4, units = "in")
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
#' @param view_sample The name of a sample to plot (as character)('FALSE' by default).
#' @param output_dir The name of a directory to write all plots to. Will create one if nonexistent.
#' @param n_threads Number of threads to use ('1' by default).
#' @param add_bins Add points for bins ('FALSE' by default). If "ref" the bins normalized copy number (divided by average (median by default)). If "alt" the bins overlapping with anchors are also added at the starting position of the reference anchor. If "FALSE" no bins are added.
#' @param prog_ploidy Ploidy of the progenitors (Assumed to be equal. '2' by default).
#' @param color_map A vector of colors for each progenitor. If "FALSE" the colors are choosen using rainbow().
#' @param specific_chr A vector of characters indicating which chromosomes to plot (plots all by default).
#' @param return_list Logical: return a list of plots if view_sample.
#'
#' @return Either nothing or a list of plots.
#' @export
#'
plot_alignment <- function(heal_list, alignment, view_sample = FALSE, output_dir = FALSE, n_threads = 1, add_bins = FALSE, prog_ploidy = 2, color_map = FALSE, specific_chr = FALSE, return_list = FALSE) {
  
  if (!add_bins %in% c(FALSE, "ref", "alt")) {
    cat("ERROR: Please input a valid 'add_bins' value. Allowed are: FALSE, 'ref' and 'alt'. Exiting..")
    return()
  }

  polyploid_samples <- names(alignment)

  sample_averages <- get_sample_stats(heal_list)

  progenitors <- names(heal_list)

  if (isFALSE(color_map)) {
    color_map <- grDevices::rainbow(length(progenitors), s = 0.7)
  } else if (length(color_map) != length(progenitors)) {
    cat("ERROR: Custom color_map is not of correct length. It should match the number of subgenomes. Exiting..")
    return()
  }
  names(color_map) <- progenitors

  if (view_sample != FALSE) {
    if (sum(polyploid_samples == view_sample) == 0) {
      cat("ERROR: Sample name not recognized (or not polyploid) for viewing of alignment. Exiting..")
      return()
    }

    if (output_dir == FALSE) {
      cat(paste0("Plotting for ", view_sample, ". \n"))
    } else {
      cat(paste0("Saving ", view_sample, "to ", output_dir, "."))
    }

    polyploid_samples <- view_sample
  } else if (output_dir == FALSE) {
    cat("ERROR: no output directory and no 'view_sample' set. One must be set. Exiting..")
    return()
  } else {
    cat(paste0("Plotting all samples and chromosomes to ", output_dir, "."))
  }

  foreach::foreach(smp = polyploid_samples) %do% {
    if (output_dir != FALSE) {
      smp_dir <- paste0(output_dir, "/", smp, "/")
      dir.create(smp_dir, showWarnings = FALSE, recursive = TRUE)
    }

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
          out_file <- paste0(ref_dir, "/", chr, ".png")
        }

        ref_chr_col_name <- paste0("chr_", ref)
        ref_start_col_name <- paste0("start_", ref)
        ref_end_col_name <- paste0("end_", ref)
        ref_cn_col_name <- paste0("cn_", ref)

        which_rows_aln_dt <- alignment[[smp]][[ref_chr_col_name]] == chr

        x_line <- (alignment[[smp]][[ref_start_col_name]][which_rows_aln_dt] + alignment[[smp]][[ref_end_col_name]][which_rows_aln_dt]) / 2
        x_vec_line <- x_line
        y_vec_line <- alignment[[smp]][[ref_cn_col_name]][which_rows_aln_dt]
        subgnm_group <- rep(ref, length(x_vec_line))

        for (alt in alt_gnms) {
          alt_cn_col_name <- paste0("cn_", alt)
          y_alt <- alignment[[smp]][[alt_cn_col_name]][which_rows_aln_dt]

          subgnm_group <- c(subgnm_group, rep(alt, length(y_alt)))
          x_vec_line <- c(x_vec_line, x_line) # same coordinates
          y_vec_line <- c(y_vec_line, y_alt)
        }

        lines_df <- data.frame(start = x_vec_line, copy = y_vec_line, subgenome = subgnm_group)

        if (add_bins != FALSE) {
          which_rows_bins_dt <- heal_list[[ref]]$bins$chr == chr

          x_pts <- heal_list[[ref]]$bins$start[which_rows_bins_dt]
          x_vec_pts <- x_pts
          y_vec_pts <- heal_list[[ref]]$bins[[smp]][which_rows_bins_dt]
          y_vec_pts <- (y_vec_pts / sample_averages[[smp]]) * prog_ploidy
          subgnm_group <- rep(ref, length(y_vec_pts))


          if (add_bins == "alt") {
            for (alt in alt_gnms) {
              alt_bin_col_name <- paste0("bin_index_", alt)
              ref_start_col_name <- paste0("start_", ref)

              cn_bindex_dt_list <- apply(alignment[[smp]][which_rows_aln_dt, ], 1, function(row) {
                bindex_vec <- as.numeric(unlist(strsplit(row[[alt_bin_col_name]], ",")))
                start_vec <- rep(row[[ref_start_col_name]], length(bindex_vec))
                return(data.table::data.table(bin_index = bindex_vec, ref_start = start_vec))
              })

              bin_cn_dt <- unique(data.table::rbindlist(cn_bindex_dt_list))

              merge_dt <- data.table::data.table(chr = heal_list[[alt]]$CN$chr[bin_cn_dt$bin_index], start = heal_list[[alt]]$CN$start[bin_cn_dt$bin_index])

              merge_dt <- merge(merge_dt, heal_list[[alt]]$bins, sort = FALSE)

              x_pts_alt <- as.numeric(bin_cn_dt$ref_start)
              y_pts_alt <- merge_dt[[smp]]
              y_pts_alt <- (y_pts_alt / sample_averages[[smp]]) * prog_ploidy

              x_vec_pts <- c(x_vec_pts, x_pts_alt)
              y_vec_pts <- c(y_vec_pts, y_pts_alt)
              subgnm_group <- c(subgnm_group, rep(alt, length(y_pts_alt)))
            }
          }

          pts_df <- data.frame(start = x_vec_pts, points = y_vec_pts, subgenome = subgnm_group)

          bin_plot <- ggplot2::ggplot() +
            ggplot2::geom_line(data = lines_df, ggplot2::aes(x = start, y = copy, color = subgenome), linewidth = 2) +
            ggplot2::geom_point(data = pts_df, ggplot2::aes(x = start, y = points, color = subgenome), size = 1, alpha = 0.1) +
            ggplot2::theme_minimal() +
            ggplot2::scale_color_manual(values = color_map) +
            ggplot2::ylim(0, 8) +
            ggplot2::labs(title = paste(chr, ref, smp), x = "Position", y = "Copy Number") +
            ggplot2::theme_bw()

          if (output_dir != FALSE) {
            ggplot2::ggsave(filename = out_file, bin_plot, device = "png", width = 6, height = 4, units = "in")
          } else {
            return(bin_plot)
          }
        } else {
          aln_plot <- ggplot2::ggplot() +
            ggplot2::geom_line(data = lines_df, ggplot2::aes(x = start, y = copy, color = subgenome), linewidth = 2) +
            ggplot2::theme_minimal() +
            ggplot2::scale_color_manual(values = color_map) +
            ggplot2::ylim(0, 8) +
            ggplot2::labs(title = paste(chr, smp, ref), x = "Position", y = "Copy Number") +
            ggplot2::theme_bw()

          if (output_dir != FALSE) {
            ggplot2::ggsave(filename = out_file, aln_plot, device = "png", width = 6, height = 4, units = "in")
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
#' @param color_vec A vector of colors for each copy number class from 0 to the ploidy of the progenitors times the number of progenitors ('FALSE' by default). If "FALSE" the colors are choosen using rainbow().
#' @param prog_ploidy Ploidy of the progenitors (Assumed to be equal. '2' by default)
#' @param n_threads Number of threads to use ('1' by default).
#'
#' @return Nothing. Plots are shown and/or saved to output_dir.
#' @export
#'
#' @importFrom data.table :=
plot_densities <- function(densities, view_sample = FALSE, output_dir = FALSE, show_discordant = FALSE, heal_list = FALSE, alignment = FALSE, corrected_alignment = FALSE, ylim_max = FALSE, color_vec = FALSE, prog_ploidy = 2, n_threads = 1) {
  is_align_and_count_data <- is.list(alignment) & is.list(heal_list)
  if (show_discordant == TRUE & !is_align_and_count_data) {
    cat("ERROR: show_discordant is TRUE but no valid heal_list or aligment provided. Exiting..")
    return()
  }

  polyploid_samples <- names(densities)
  progenitors <- names(heal_list)

  if (view_sample != FALSE) {
    n_threads <- 1
    if (sum(polyploid_samples == view_sample) == 0) {
      cat("ERROR: Sample name not recognized (or not polyploid) for viewing of alignment. Exiting..")
      return()
    }

    cat(paste0("Plotting for ", view_sample, ". \n"))
    polyploid_samples <- view_sample
  } else if (output_dir == FALSE) {
    cat("ERROR: no output directory and no 'view_sample' set. One must be set. Exiting..")
    return()
  } else {
    cat(paste0("Plotting all samples and chromosomes to ", output_dir, "."))
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
        cat("ERROR: No ylim_max and no heal_list. Set at least one. Exiting..")
        return()
      } else {
        counts_vec <- unlist(lapply(heal_list, function(dt) {
          return(dt$bins[[smp]])
        }))
        ylim_max <- mean(counts_vec, na.rm = TRUE) + 3 * stats::sd(counts_vec, na.rm = TRUE)
      }
    }

    cn_labels <- names(densities[[smp]])

    if (isFALSE(color_vec)) {
      color_vec <- grDevices::rainbow(n = length(cn_labels), s = 0.7, v = 1)
    } else {
      if (length(cn_labels) != length(color_vec)) {
        cat("ERROR: color vector length not matching number of copy number categories. Exiting..")
        return()
      }
    }

    names(color_vec) <- cn_labels

    graphics::contour(densities[[smp]][[cn_labels[1]]], ylim = c(-5, ylim_max), col = color_vec[cn_labels[1]], main = smp)

    for (i in 2:length(cn_labels)) {
      graphics::contour(densities[[smp]][[cn_labels[i]]], ylim = c(0, ylim_max), col = color_vec[cn_labels[i]], add = TRUE)
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

        merge_dt <- merge(heal_list[[prog]]$bins, merge_dt)

        col_keep <- c("gc_content", smp)

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
          merge_dt <- merge(heal_list[[prog]]$bins, merge_dt)
          merge_dt <- merge(heal_list[[prog]]$CN, merge_dt, by = c("chr", "start", "gc_content"), suffixes = c("_CN", "_count"))
          col_keep <- c("gc_content", paste0(smp, c("_CN", "_count")))
          merge_dt <- merge_dt[, ..col_keep]
          data.table::setnames(merge_dt, old = col_keep, c("gc_content", "cn_original", "count"))
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


#' Plot linear regression of copy number at anchors between each pair of progenitors
#'
#' @param alignment A heal alignment object created with get_heal_alignment().
#' @param view_samples A vector of sample names to plot (as character)('FALSE' by default).
#' @param output_dir The name of a directory to write all plots to. Will create one if nonexistent.
#' @param color The color of the points ("blue4" by default).
#' @param alpha The transparency of the points (0.1 by default).
#' @param width The jitter of points along the x axis.
#' @param height The jitter of points along the x axis.
#'
#' @return Nothing. Plots are shown and/or saved to output_dir.
#' @export
#'
plot_linear_alignment <- function(alignment, view_samples = FALSE, output_dir = FALSE, color = "blue4", alpha = 0.1, width = 0.2, height = 0.2) {
  polyploid_samples <- names(alignment)


  if (!isFALSE(view_samples)) {
    if (length(intersect(polyploid_samples, view_samples)) == 0) {
      cat("ERROR: Sample names not recognized (or not polyploid) for viewing of alignment. Exiting..")
      return()
    }

    cat(paste("Plotting for:", view_samples, " \n"))
    polyploid_samples <- view_samples
  } else if (isFALSE(output_dir)) {
    cat("ERROR: no output directory and no 'view_samples' set. One must be set. Exiting..")
    return()
  } else {
    cat(paste0("Plotting all samples to ", output_dir, "."))
  }


  for (smp in polyploid_samples) {
    cn_cols <- grep("cn_", colnames(alignment[[smp]]), value = TRUE)
    all_pairs <- utils::combn(cn_cols, 2)

    for (i in 1:ncol(all_pairs)) {
      if (output_dir != FALSE) {
        dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
        out_file <- paste0(output_dir, "/", pair[1], "_vs_", pair[2], "_linear.png")
      }

      pair <- all_pairs[, i]
      progenitors <- sub("cn_", "", pair)

      input_dt <- data.table::data.table(x = alignment[[smp]][[pair[1]]], y = alignment[[smp]][[pair[2]]])

      plot <- ggplot2::ggplot(input_dt, ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_jitter(color = color, alpha = alpha, width = width, height = height) +
        ggplot2::labs(
          title = smp,
          x = paste("Infered CN in", progenitors[1]), y = paste("Infered CN in", progenitors[2])
        ) +
        ggplot2::theme_minimal()

      if (output_dir != FALSE) {
        ggplot2::ggsave(filename = out_file, plot, device = "png", width = 6, height = 4, units = "in")
      } else {
        print(plot)
      }
    }
  }
}


#' Plot pairwise heatmap
#'
#' @param alignment A heal alignment object created with get_heal_alignment().
#' @param view_samples A vector of sample names to plot (as character)('FALSE' by default).
#' @param output_dir The name of a directory to write all plots to. Will create one if nonexistent.
#'
#' @return Nothing. Plots are shown and/or saved to output_dir.
#' @export
#'
plot_heal_heat_map <- function(alignment, view_samples = FALSE, output_dir = FALSE) {
  polyploid_samples <- names(alignment)


  if (!isFALSE(view_samples)) {
    if (length(intersect(polyploid_samples, view_samples)) == 0) {
      cat("ERROR: Sample names not recognized (or not polyploid) for viewing of alignment. Exiting..")
      return()
    }

    cat(paste("Plotting for:", view_samples, " \n"))
    polyploid_samples <- view_samples
  } else if (isFALSE(output_dir)) {
    cat("ERROR: no output directory and no 'view_samples' set. One must be set. Exiting..")
    return()
  } else {
    cat(paste0("Plotting all samples to ", output_dir, "."))
  }


  for (smp in polyploid_samples) {
    cn_cols <- grep("cn_", colnames(alignment[[smp]]), value = TRUE)
    all_pairs <- utils::combn(cn_cols, 2)

    for (i in 1:ncol(all_pairs)) {
      if (output_dir != FALSE) {
        dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
        out_file <- paste0(output_dir, "/", pair[1], "_vs_", pair[2], "_heat.png")
      }

      pair <- all_pairs[, i]
      progenitors <- sub("cn_", "", pair)

      input_dt <- data.table::data.table(x = alignment[[smp]][[pair[1]]], y = alignment[[smp]][[pair[2]]])

      count_table <- as.data.frame(table(input_dt))
      colnames(count_table) <- c("x", "y", "count")
      count_table$x <- as.factor(count_table$x)
      count_table$y <- as.factor(count_table$y)

      plot <- ggplot2::ggplot(count_table, ggplot2::aes(x = x, y = y, fill = count)) +
        ggplot2::geom_tile(color = "black") +
        ggplot2::scale_fill_gradient(low = "blue", high = "red", na.value = "grey50") +
        ggplot2::labs(
          title = paste0("Heatmap of Counts for Homoeolog Pairs in ", smp),
          x = paste("Infered CN in", progenitors[1]), y = paste("Infered CN in", progenitors[2])
        ) +
        ggplot2::theme_minimal()

      if (output_dir != FALSE) {
        ggplot2::ggsave(filename = out_file, plot, device = "png", width = 6, height = 4, units = "in")
      } else {
        print(plot)
      }
    }
  }
}

utils::globalVariables(c("a", "..col_keep", "cn", "i", "status", "..bin_index_col", "..which_relevant", "..cn_col", "bin_start", "anchor_start", "anchor_end", "gene_id", "ref", "copy", "subgenome", "points", "progenitor", "counts", "count", "is_allowed", "x", "y"))
