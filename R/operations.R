compare_heal_list <- function(heal_list_1, heal_list_2){
  
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
    }
  }else{
    stop("Input healr lists don't have the same progenitor sets. Exiting..")
  }
  
  get_sample_stats(heal_list_1, )
}