#' Title
#'
#' @param x expression matrix as data.table (recommended) or data.frame object whose first column (not the rownames!) holds feature names
#' @param y ditto.
#' @param z ditto.
#' @param output.dir
#' @param max.iter
#' @param max.K
#' @param rep.runs
#' @param pItem
#' @param pFeature
#' @param clusterAlg
#' @param distance
#' @param cc.seed
#'
#' @return
#' @export
#'
#' @examples
BIConsensusCluster <- function(x, y, z, output.dir, max.iter = 20, max.K = 7, rep.runs = 1000, pItem=0.8, pFeature=1, clusterAlg="hc", distance="euclidean", cc.seed=5000){
  invisible(sapply(list(x, y, z), function(x) setnames(x, 1, 'ID_REF')))
  iteration <- 1
  dir.create(output.dir)
  setwd(output.dir)
  gene.sig <- x[,1]

  while(iteration <= max.iter){
    print(paste(date(), iteration, sep=" -- start iteration: "))
    platforms <- c('x', 'y', 'z')
    run.dir <- vapply(platforms,
                      function(x) paste(x, as.character(iteration), sep = ""),
                      "Yu Fat is handsome")
    x.cc <- ConsensusClusterPlus(as.matrix(x[gene.sig, on = 'ID_REF'][,-1]),maxK=max.K,reps=rep.runs,pItem=pItem,pFeature=pFeature, title=run.dir['x'] ,clusterAlg=clusterAlg,distance=distance,seed=cc.seed,plot="png")
    y.cc <- ConsensusClusterPlus(as.matrix(y[gene.sig, on = 'ID_REF'][,-1]),maxK=max.K,reps=rep.runs,pItem=pItem,pFeature=pFeature, title=run.dir['y'],clusterAlg=clusterAlg,distance=distance,seed=cc.seed,plot="png")
    z.cc <- ConsensusClusterPlus(as.matrix(z[gene.sig, on = 'ID_REF'][,-1]),maxK=max.K,reps=rep.runs,pItem=pItem,pFeature=pFeature, title=run.dir['z'],clusterAlg=clusterAlg,distance=distance,seed=cc.seed,plot="png")

    iteration<- iteration + 1
  }
}
