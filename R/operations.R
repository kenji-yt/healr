#' Find the difference in copy number between shared samples of two heal lists. 
#'
#' @param heal_list_1 List in heal format (such as output from count_heal_data()).
#' @param heal_list_2 List in heal format (such as output from count_heal_data()).
#' @param n_threads Number of threads to use ('1' by default). It parallelize for each shared sample. 
#' The maximum number used is therefore the number of shared samples.
#' @param abs Return the absolute value of the differences ('FALSE' by default).
#'
#' @returns A heal list with one data table per subgenome holding the difference between shared samples (heal_list_1 - heal_list_2).
#' @export
#'
#' @examples
compare_heal_list <- function(heal_list_1, heal_list_2, n_threads, abs = FALSE){
  
  var_name_1 <- substitute(heal_list_1)
  var_name_2 <- substitute(heal_list_2)
  
  cn_is_null_1 <- is.null(unlist(lapply(heal_list_1, function(list) {
    list$CN
  })))
  cn_is_null_2 <- is.null(unlist(lapply(heal_list_2, function(list) {
    list$CN
  })))
  
  if(cn_is_null_1 | cn_is_null_2){
    stop("At least one list has no CN data. Exiting..")
  }
  
  progs_1 <- names(heal_list_1)
  progs_2 <- names(heal_list_2)
  
  if(length(progs_1)==length(progs_2)){
    diffs <- sum(progs_1!=progs_2)
    if(diffs!=0){
      stop("Input healr lists don't have the same progenitor sets. Exiting..")
    }else{
      progs <- progs_1
    }
  }else{
    stop("Input healr lists don't have the same progenitor sets. Exiting..")
  }
  
  samples_1 <- names(get_sample_stats(heal_list_1))
  samples_2 <- names(get_sample_stats(heal_list_2))
  
  common_samples <- intersect(samples_1, samples_2)
  if(length(common_samples)==0){
    stop("No common samples to compare. Exiting..")
  }
  
  for(prog in progs){
    if(nrow(heal_list_1[[prog]]$CN)!=nrow(heal_list_2[[prog]]$CN)){
      stop(paste0("Unequal number of rows in the CN table of ", prog, ". Exiting.."))
    }
  }
  
  diff_list <- foreach::foreach(prog = progs)%do%{
    
    general_info_dt <- heal_list_1[[prog]]$CN[,1:4]
    
    doParallel::registerDoParallel(n_threads)
    diff_per_smp <- foreach::foreach(smp = common_samples)%dopar%{
      diff_vec <- heal_list_1[[prog]]$CN[[smp]] - heal_list_2[[prog]]$CN[[smp]]
      if(abs==TRUE){
        diff_vec <- abs(diff_vec)
      }
      return(diff_vec)
    }
    doParallel::stopImplicitCluster()
    names(diff_per_smp) <- common_samples
    
    diff_dt_prog <- data.table::as.data.table(diff_per_smp)
    
    diff_dt <- cbind(general_info_dt, diff_dt_prog)
    return(list(diff=diff_dt))
  }
  names(diff_list) <- progs
  
  return(diff_list)
}
