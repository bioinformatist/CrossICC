#' Title
#'
#' @param consencus.result
#' @param cutoff
#' @param maxK
#'
#' @return
#' @export
#'
#' @examples
derive.clusternum <- function(consencus.result, cutoff = 0.05, maxK){
  for(k in 2:maxK){
    A=cal.auc(consencus.result,k)
    if(k == 2){
      deltA=A
      preA=A
    } else deltA=(A-preA)/preA
    if(deltA<=cutoff) return(k-1)
    preA=A
  }
}
