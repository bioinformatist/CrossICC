#' Title
#'
#' @param x matrix (data.frame) with sample names as column names and feature names as row names.
#' @param y ditto.
#' @param z ditto.
#' @param output.dir the results' output directory.
#' @param max.iter the maximum number of iterations.
#' @param max.K the maximum number of cluster center.
#' @param rep.runs number of subsamples during clustering.
#' @param pItem proportion of items to sample during clustering.
#' @param pFeature proportion of features to sample during clustering.
#' @param clusterAlg cluster algorithm. Could be 'hc' heirarchical (hclust),
#' 'pam' for paritioning around medoids, 'km' for k-means upon data matrix,
#' 'kmdist' for k-means upon distance matrices (former km option),
#' or a function that returns a clustering.
#' @param distance Could be 'pearson': (1 - Pearson correlation),
#' 'spearman' (1 - Spearman correlation),
#' 'euclidean', 'binary', 'maximum', 'canberra', 'minkowski" or custom distance function.
#' @param cc.seed sets random seed for reproducible results.
#' @param filter.cutoff cutoff value during mad filtering.
#' @param fdr.cutoff cutoff value during fdr filtering.
#' @param cluster.cutoff cutoff value during determining cluster numbers.
#' @param ebayes.cutoff p-value cutoff when select differentially expressed probes.
#'
#' @return A nested list with iteration time as its name and list containing consensus cluster,
#' gene signature and balanced cluster as its value.
#' @export
#'
#' @seealso \code{\link[ConsensusClusterPlus]{ConsensusClusterPlus}}
#'
#' @examples
#' BIConsensusCluster(x, y, z, output.dir = 'handsome_Yu_Fat', max.iter = 1)
BIConsensusCluster <- function(x, y, z, filter.cutoff = 0.5, fdr.cutoff = 0.1, output.dir, max.iter = 20, max.K = 6, rep.runs = 1000,
                               pItem=0.8, pFeature=1, clusterAlg="hc", distance="euclidean",
                               cc.seed=5000, cluster.cutoff = 0.05, ebayes.cutoff = 0.1){
  dir.create(output.dir)
  setwd(output.dir)

  # Merge, filter and scale here
  cat(paste(date(), '--', 'Pre-processing data'), '\n')
  mfs.list <- m.f.s(x, y, z, fdr.cutoff = fdr.cutoff, filter.cutoff = filter.cutoff)
  iteration <- 1
  platforms <- list(x = mfs.list[[2]][[1]], y = mfs.list[[2]][[2]], z = mfs.list[[2]][[3]])
  gene.sig <- rownames(platforms$x)

  while(iteration <= max.iter){
    cat(paste(date(), iteration, sep=" -- start iteration: "), '\n')
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

    xyz.sig <- list(x = platforms$x[gene.sig,],
                    y = platforms$y[gene.sig,],
                    z = platforms$z[gene.sig,])
    balanced.cluster <- balance.cluster(xyz.sig,
                                        cc = cc, cluster.cutoff = cluster.cutoff,
                                        max.K = max.K, plot = TRUE, iter = iteration)

    gene.sig.all <- lapply(names(xyz.sig),
                           function(x) rownames(ebayes(xyz.sig[[x]],
                                                       balanced.cluster[[2]][[x]],
                                                       cutoff = ebayes.cutoff)[[1]]))
    pre.gene.sig <- gene.sig
    gene.sig <- com.feature(unlist(gene.sig.all), method = 'merge')

    if(!isTRUE(all.equal(pre.gene.sig, gene.sig)) && isTRUE(all.equal(sort(pre.gene.sig), sort(gene.sig)))){
      break
    }

    result[[iteration]] <- list(consensus.cluster = cc,
                                gene.signature = gene.sig,
                                balanced.cluster = balanced.cluster)

    iteration<- iteration + 1
  }
  result
}
