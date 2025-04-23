#' Write a heal list to a directory
#'
#' @param heal_list List in heal format (such as output from count_heal_data()).
#' @param output_dir The name of a directory to write the heal list to. Fails if the directory exists.
#' @param ... Any argument you wish to pass to fwrite (see data.table::fwrite documentation). Writes as csv by default. /!\ Warning: Changing from default might cause issues in read_heal_list() function.
#'
#' @returns Nothing. Writes the data into a directory. 
#' @export
#'
write_heal_list <- function(heal_list, output_dir, ...){
  
  if(dir.exists(output_dir)) {
    stop(paste("Output directory",output_dir,"already exists!"))
    
  } else {
    
    progenitors <- names(heal_list)
    
    for(prog in subgenomes){
      
      prog_dir <- paste0(output_dir,"/",prog)
      
      dir.create(prog_dir, recursive = TRUE)
      
      dt_types <- names(heal_list[[prog]])
      
      for(dt_t in dt_types){
        
        dt <- heal_list[[prog]][[dt_t]]
        
        file_path <- paste0(prog_dir, "/", dt_t,".dt")
        data.table::fwrite(dt, file=file_path, ...)
        
      }
      cat(paste0("Done for ",prog,".\n"))
    }
  }
}


#' Read a heal list from a directory
#' 
#' @param heal_list_directory A directory structured after a heal list (created by write_heal_list()).
#' @param ... Any argument you which to pass to fread (see data.table::fread documentation).
#'
#' @returns A heal list object i.e. a list with one element per progenitor containing at least a data table with counts in bins for each sample and GC and mappability for each bin.
#' @export
#'
#' @examples
read_heal_list <- function(heal_list_directory, ...){
  
  out_heal_list <- list()
  
  progenitors <- list.dirs(path = heal_list_directory, full.names = FALSE, recursive = FALSE)
  progenitor_dirs <- paste0(heal_list_directory,"/",progenitors)
  names(progenitor_dirs) <- progenitors
  
  for(prog in progenitors){
    
    prog_list <- list()
    
    file_names <- list.files(progenitor_dirs[prog], full.names = FALSE)
    file_paths <- paste0(progenitor_dirs[prog],"/",file_names)
    dt_names <- gsub(".dt","", file_names)
    names(file_paths) <- dt_names
    
    for(dt_name in dt_names){
      prog_list[[dt_name]] <- data.table::fread(file=file_paths[dt_name], ...)
    }
    
    out_heal_list[[prog]] <- prog_list
    
  }
  
  return(out_heal_list)
  
}


#' Write a heal copy number summary to a directory
#'
#' @param cn_summary A list object created with summarize_cn() 
#' @param output_dir The name of a directory to write the summary list to. Fails if the directory exists.
#' @param ... Any argument you wish to pass to fwrite (see data.table::fwrite documentation). Writes as csv by default. /!\ Warning: Changing from default might cause issues in read_cn_summary() function.
#'
#' @returns Nothing. Writes the data into a directory. 
#' @export
#'
#' @examples
write_cn_summary <- function(cn_summary, output_dir, ...){
  
  if(dir.exists(output_dir)) {
    stop(paste("Output directory",output_dir,"already exists!"))
    
  } else {
    
    dir.create(output_dir, recursive = TRUE)
    
    samples <- names(cn_summary)
    for(smp in samples){
    
      smp_dir <- paste0(output_dir,"/",smp)
      dir.create(smp_dir)
      
      file_path_total <- paste0(smp_dir, "/total_count_table.dt")
      data.table::fwrite(cn_summary[[smp]]$total_count_table, file = file_path_total, ...)
      
      progenitors <- setdiff(names(cn_summary[[smp]]),"total_count_table")
      
      for(prog in progenitors){
        
        smp_prog_dir <- paste0(smp_dir,"/",prog)
        dir.create(smp_prog_dir)
        
        dt_name <- paste0("total_count_table_", prog)
        filename <- paste0(dt_name, ".dt")
        file_path_sub <- paste0(smp_prog_dir,"/",filename)
        data.table::fwrite(cn_summary[[smp]][[prog]][[dt_name]], file = file_path_sub, ...)
        
        chromosomes <- setdiff(names(cn_summary[[smp]][[prog]]), dt_name)
        
        for(chromo in chromosomes){
          
          chromo_dir <- paste0(smp_prog_dir,"/",chromo)
          dir.create(chromo_dir)
          
          count_matrix <- cn_summary[[smp]][[prog]][[chromo]]$count_table
          count_table <- data.table::data.table(copy_number=colnames(count_matrix), counts=count_matrix[1,], percentages=count_matrix[2,])
          path_count_dt <- paste0(chromo_dir, "/count_table.dt")
          data.table::fwrite(count_table, file = path_count_dt, ...)
          
          rle <- cn_summary[[smp]][[prog]][[chromo]]$run_length_encoding 
          rle_dt <- data.table::data.table(lengths=rle$lengths, values=rle$values)
          path_rle_dt <- paste0(chromo_dir, "/run_length_encoding.dt")
          data.table::fwrite(rle_dt, file = path_rle_dt, ...)
          
        }
      }
    }
    cat(paste0("Done for ", smp, ".\n"))
  }
}


#' Read a heal copy number summary from a directory
#'
#' @param cn_summary_directory A directory structured after a CN summary (created by write_cn_summary()).
#' @param ... Any argument you which to pass to fread (see data.table::fread documentation).
#'
#' @returns A CN summary i.e. a list object with information about copy number for each sample. Same as the output of summarize_cn().
#' @export
#'
#' @examples
read_cn_summary <- function(cn_summary_directory, ...){
  
  out_cn_summary <- list()
  
  samples <- list.dirs(path = cn_summary_directory, full.names = FALSE, recursive = FALSE)
  sample_dirs <- paste0(cn_summary_directory,"/",samples)
  names(sample_dirs) <- samples
  
  for(smp in samples){
    
    smp_list <- list()
    
    file_names <- list.files(sample_dirs[smp], full.names = FALSE)
    
    smp_list$total_count_table <- data.table::fread(file=paste0(sample_dirs[smp], "/total_count_table.dt"), header=TRUE, ...)
    
    progenitors <- setdiff(file_names, "total_count_table.dt")
    smp_paths <- paste0(sample_dirs[smp], "/", progenitors)
    names(smp_paths) <- progenitors
    
    for(prog in progenitors){
      
      prog_list <- list()
      
      prog_file_names <- list.files(smp_paths[prog], full.names = FALSE)
      
      total_dt_name <- paste0("total_count_table_", prog, ".dt")
      
      prog_list[[total_dt_name]] <- data.table::fread(file=paste0(smp_paths[prog], "/", total_dt_name), header=TRUE, ...)
      
      chromosomes <- setdiff(prog_file_names, total_dt_name)
      
      chr_paths <- paste0(smp_paths[prog], "/", chromosomes)
      names(chr_paths) <- chromosomes
      
      for(chromo in chromosomes){
        
        chromo_list <- list()
        
        chromo_list$count_table <- data.table::fread(file=paste0(chr_paths[chromo], "/count_table.dt"), header=TRUE, ...)
        
        rle_dt <- data.table::fread(file = paste0(chr_paths[chromo],"/run_length_encoding.dt"), header=TRUE, ...)
        run_length_encoding <- list(values=rle_dt$values, lengths=rle_dt$lengths)
        class(run_length_encoding) <- "rle"
        
        chromo_list$run_length_encoding <- run_length_encoding
        
        prog_list[[chromo]] <- chromo_list
      }
      
      smp_list[[prog]] <- prog_list
    }
    
    out_cn_summary[[smp]] <- smp_list
    
  }
  
  return(out_cn_summary)
  
}

write_heal_alignment <- function(alignment, output_dir, ...){
  
  if(dir.exists(output_dir)) {
    stop(paste("Output directory", output_dir, "already exists!"))
    
  } else {
    
    dir.create(output_dir, recursive = TRUE)
    
    samples <- names(alignment)
    for(smp in samples){
      
      smp_dt_path <- paste0(output_dir,"/",smp,".dt")
      data.table::fwrite(alignment[[smp]], file = smp_dt_path, ...)
      
    }
    cat(paste0("Done for ", smp, ".\n"))
  }
}

read_heal_alignment <- function(alignment_directory, ...){
  
  out_alignment <- list()
  
  sample_files <- list.files(path = alignment_directory, full.names = FALSE, recursive = FALSE)
  sample_names <- gsub(".dt", "", sample_files)
  names(sample_files) <- sample_names
  
  for(smp in sample_names){
    
    out_alignment[[smp]] <- data.table::fread(sample_files[smp], header=TRUE, ...)
  }
  return(out_alignment)
}

write_alignment_summary <- function(alignment_summary, output_dir, formating){
  
  if (dir.exists(dir_path)) {
    stop("Output directory already exists!")
    
  } else {
    dir.create(output_dir, recursive = TRUE)
    
    subgenomes <- names(heal_list)
  }
}