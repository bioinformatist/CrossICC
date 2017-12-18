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
      ΔA=A
      pre.A=A
    } else ΔA=(A-pre.A)/pre.A
    if(ΔA<=cutoff) return(k-1)
    pre.A=A
  }
}
