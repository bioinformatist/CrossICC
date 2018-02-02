derive.clusternum <- function(consencus.result, cutoff = 0.05, maxK){
  for(k in 2:maxK){
    A <- cal.auc(consencus.result, k)
    if(k == 2){
      delta.A <- A
      pre.A <- A
    } else delta.A <- (A-pre.A) / pre.A
    if (delta.A <= cutoff) {
      # Return value must be larger than 2!
      if ((k - 1) < 2) {
        return(2)
      } else {
        return(k - 1)
      }
    }
    pre.A <- A
  }
}
