#' Title
#'
#' @param x matrix (data.frame) with sample names as column names and feature names as row names.
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
#' @param filter.cutoff
#' @param fdr.cutoff
#'
#' @return
#' @export
#'
#' @examples
BIConsensusCluster <- function(x, y, z, filter.cutoff = 0.5, fdr.cutoff = 0.001, output.dir, max.iter = 20, max.K = 6, rep.runs = 1000,
                               pItem=0.8, pFeature=1, clusterAlg="hc", distance="euclidean",
                               cc.seed=5000){
  dir.create(output.dir)
  setwd(output.dir)
  mfs.list <- m.f.s(x, y, z)
  iteration <- 1
  gene.sig <- rownames(mfs.list[[2]][[1]])

  # Merge, filter and scale here
  platforms <- list(x = mfs.list[[2]][[1]], y = mfs.list[[2]][[2]], z = mfs.list[[2]][[3]])

  while(iteration <= max.iter){
    print(paste(date(), iteration, sep=" -- start iteration: "))
    # Here a named vector is needed (for dir names), so (v)applys is necessary
    run.dir <- vapply(names(platforms),
                      function(x) paste(x, iteration, sep = ""),
                      "Yu Fat is handsome")
    cc <- vapply(names(platforms),
                 function(x) list(x = ConsensusClusterPlus::ConsensusClusterPlus(platforms[[x]][gene.sig,],
                                                                                 maxK=max.K, reps=rep.runs, pItem=pItem,
                                                                                 pFeature=pFeature, title=run.dir[x],
                                                                                 clusterAlg=clusterAlg, distance=distance,
                                                                                 seed=cc.seed, plot="png")),
                 # ConsensusClusterPlus returns a list of 7 elements, but we need a nested list
                 list(x = rep(list('fuck'), 7)))
    iteration<- iteration + 1
  }
  return(cc)
}
