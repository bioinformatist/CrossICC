#' CrossICC: Automatically Aggregating and Summarizing Bioinformatics Results for Interactive Report.
#'
#' The CrossICC package provides three categories of important functions:
#' foo, bar and baz.
#'
#' @section CrossICC functions:
#' The foo functions ...
#'
#' @docType package
#' @name CrossICC
#' @import data.table MASS MergeMaid grDevices stats
#' @importFrom methods new
#' @importFrom utils head
NULL

#' Title The Main Function of the package
#'
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
#' @param study.names a vector containing all study names
#' @param ... all datasets (matrices is better)
#' @param method 'finer' (based on correlation between sample expression values and centroids) or 'balanced' (based on correlation of centroids).
#' Which super-cluster strategy should be used?
#' @param ebayes.mode 'up' or 'both'. Choose only up-regulated genes or all differentially expressed genes.
#' @param use.shiny If TRUE, a shiny app will appear after running this main function.
#'
#' @return A nested list with iteration time as its name and list containing consensus cluster,
#' gene signature and balanced cluster as its value.
#'
#' @export
#'
#' @seealso \code{\link[ConsensusClusterPlus]{ConsensusClusterPlus}}
#'
#' @examples
#' \donttest{
#' # It takes too long time for running code below, so ignore them in R CMD check.
#' CrossICC(example.matrices, max.iter = 20, use.shiny = FALSE)
#' CrossICC(example.matrices, output.dir = 'handsome_Yu_Fat', max.iter = 20)
#' fuck <- CrossICC(datalist, max.iter = 5, use.shiny = FALSE, method = "balanced",fdr.cutoff = 0.5, ebayes.cutoff = 0.01)
#' }
CrossICC <- function(..., study.names, filter.cutoff = 0.5, fdr.cutoff = 0.1, output.dir = NULL, max.iter = 20, rep.runs = 1000,
                               pItem = 0.8, pFeature = 1, clusterAlg = "hc", distance = "euclidean",
                               cc.seed = 5000, cluster.cutoff = 0.05, ebayes.cutoff = 0.1, ebayes.mode = 'up', method = 'finer', use.shiny = TRUE){
  if (max.iter < 2) warning('Result from less than 2 times iteration may not make sense at all!')

  # Check method (sometimes spelling mistake here)
  method <- match.arg(method, c("finer","balanced"))

  # Get arguments as data.table
  arg.list <- unlist(as.list(match.call())[-1])

  graphics.off()

  if(is.null(output.dir)) {
    plot.suffix = NULL
  } else {
    plot.suffix = 'png'
    dir.create(output.dir)
    setwd(output.dir)
  }

  arg <- list(...)
  result <- list()

  # It is not allowed just a matrix as input, for MergeMaid will raise an error:
  # Error in .dens.mergeExpressionSet(x = x, method = method, ...) :
  #   Number of studies in the mergeExpressionSet should not less than 2.
  # Already can auto-skip this filtering when there's only one sample.
  # if ((length(arg) == 1) & (is.element(class(arg[[1]]),"matrix"))) {
  #   stop("Number of studies should not less than 2.")
  if ((length(arg) == 1) & (is.element(class(arg[[1]]),"list"))) {
    platforms.list <- unlist(arg, recursive = FALSE)
  } else {
    platforms.list <- lapply(arg, check.eSet)
  }

  # Merge, filter and scale here
  cat(paste(date(), '--', 'Pre-processing data'), '\n')
  mfs.list <- m.f.s(platforms.list, fdr.cutoff = fdr.cutoff, filter.cutoff = filter.cutoff)
  iteration <- 1
  platforms <- mfs.list[[2]]

  if (missing(study.names)) sn <- c()

  # If study.names is not defined, use automatically generated ones instead
  if ((!is.element(class(sn), "vector")) | (!length(sn) == length(mfs.list[[2]]))) {
    study.names <- sapply(1:length(platforms), function(x) paste0('Matrix.', x))
    names(platforms) <- study.names
  } else {
    names(platforms) <- study.names
  }

  # Use the first datasets' feature names as default
  gene.sig <- rownames(platforms[[1]])

  # Use matrix save iteration time with signature number
  iter.sig <- data.table(Iteration = numeric(), Signatures = numeric())

  while(iteration <= max.iter){
    cat(paste(date(), iteration, sep=" -- start iteration: "), '\n')
    # Here a named vector is needed (for dir names), so (v)applys is necessary
    run.dir <- vapply(names(platforms),
                      function(x) paste(x, iteration, sep = "."),
                      "Yu Fat is handsome")

    pdf(NULL)
    par(mar=c(1,1,1,1))
    dev.control('enable') # enable display list
    cc <- vapply(names(platforms),
                 function(x) list(suppressMessages(ConsensusClusterPlus::ConsensusClusterPlus(platforms[[x]][gene.sig,],
                                                                             maxK = 4, reps=rep.runs, pItem=pItem,
                                                                             pFeature=pFeature, title=run.dir[x],
                                                                             clusterAlg=clusterAlg, distance=distance,
                                                                             seed=cc.seed, plot = plot.suffix))),
                 # ConsensusClusterPlus returns a list of 7 elements, but we need a nested list
                 list(rep(list('fuck'), 7)))
    dev.off()

    all.sig <- lapply(names(platforms), function(x) platforms[[x]][gene.sig,])
    names(all.sig) <- names(platforms)

    balanced.cluster <- balance.cluster(all.sig,
                                        cc = cc, cluster.cutoff = cluster.cutoff,
                                        max.K = NULL, method = method)

    ebayes.result <- lapply(names(all.sig),
                            function(x) ebayes(all.sig[[x]],
                                               switch(method,
                                                      'balanced' = balanced.cluster[[1]][[x]],
                                                      'finer' = balanced.cluster[[1]]),
                                               cutoff = ebayes.cutoff,
                                               mode = ebayes.mode))

    # Remove possible NULL element produced during ebayes step
    ebayes.result <- Filter(Negate(is.null), ebayes.result)

    gene.sig.all <- lapply(ebayes.result, function(x) rownames(x$full.m))
    geneset2gene <- lapply(ebayes.result, function(x) x$geneset2gene)

    unioned.genesets <- as.matrix(unique(data.table(do.call(rbind, geneset2gene))))

    pre.gene.sig <- gene.sig
    gene.sig <- com.feature(unlist(gene.sig.all), method = 'merge')

    merged.matrix <- do.call(cbind, all.sig)
    merged.cluster <- switch (method,
      'balanced' = do.call(c, balanced.cluster[[1]]),
      'finer' = balanced.cluster[[1]]
    )

    design <- model.matrix(~ 0 + factor(merged.cluster))
    k <- length(unique(merged.cluster))
    colnames(design) <- paste("K", 1:k, sep = "")
    K <- colnames(design)
    fit <- limma::lmFit(merged.matrix, design)

    x <- c()
    for(i in 1:(k - 1)){
      for(j in (i + 1):k){
        x <- c(x, paste(K[j], K[i], sep = "-"))
      }
    }
    for(i in 1:k){
      x <- c(x, paste(K[i], paste("(", paste(K[-i], collapse="+"), ")", "/", k-1, sep=""), sep="-"))
    }

    contrast.matrix <- limma::makeContrasts(contrasts = x,levels = design)
    fit2 <- limma::contrasts.fit(fit,contrast.matrix)
    fit2 <- limma::eBayes(fit2)
    FC.matrix <- limma::topTable(fit2, number = 20000, adjust.method = 'BH')
    sorted.genes <- rownames(FC.matrix[order(-FC.matrix[,1]),])
    gene.order <- sorted.genes[sorted.genes %in% gene.sig]

    # Old version of sorting
    # sorted.gene.list <- lapply(ebayes.result, function(x) rownames(x$full.m[rownames(x$full.m) %in% gene.sig,][order(-x$full.m[rownames(x$full.m) %in% gene.sig,][,1]),]))

    sig.num <- length(gene.sig)
    iter.sig <- rbind(iter.sig, data.table(Iteration = iteration, Signatures = sig.num))
    cat(sig.num, "genes were engaged in this iteration.\n")

    if(sig.num == 0) {
      warning("No signature was kept during last time iteration. You may change ebayes.mode to \"both\" then try again.")
      break
    }

    # Confirm those two atomic vectors are exactly equal
    # Below is for perfect match (with the same order)
    # if(isTRUE(all.equal(pre.gene.sig, gene.sig)) && isTRUE(all.equal(sort(pre.gene.sig), sort(gene.sig)))){
    if(isTRUE(all.equal(sort(pre.gene.sig), sort(gene.sig)))){
      if (iteration ==1) {
        warning("All reserved signatures after processing (merging, filtering and scaling) already reach convergence.")
      }
      break
    }

    if (iteration == max.iter) {
      warning("Max iteration time reached! May still not reach convergence.")
      break
    }
    iteration<- iteration + 1
  }

  # Just report results of last iteration (To reduce the size of RDS used by shiny)
  result <- list(# consensus.cluster = cc,  # For test only
    # all.sig = all.sig,  # For test only
    # er = ebayes.result,  # For test only
    arg.list = arg.list,  # named vector object of arguments
    platforms = platforms,  # For heatmap in shiny
    gene.signature = gene.sig,
    iter.sig = iter.sig,
    gene.order = gene.order,  # Sorted gene names by Fold-Change value, for heatmaps use
    # gene.sig.all = gene.sig.all,  # For test only
    # MDEG = gene.sig.all,  # For test only
    clusters = balanced.cluster,
    geneset2gene = geneset2gene,
    unioned.genesets = unioned.genesets)

  saveRDS(result, file = path.expand('~/CrossICC.object.rds'))
  cat("A CrossICC.object.rds file will be generated in home directory by default.
      Note that the previous file will be overridden.\n")

  cat(paste(date(), iteration - 1, sep=" -- Iteration finished! Iteration time for reaching convergence/limit: "), '\n')

  if (use.shiny) {
    warning('Result list would not be returned when use.shiny=TRUE.')
    run.shiny()
  }
  # list(all.sig = all.sig, balanced.cluster = balanced.cluster)
  result
  # all.sig
  # cc
  # list(all.sig, cc)
  # platforms
  # balanced.cluster
  # ebayes.result
}

run.shiny<-function(){
  shiny::runApp(system.file("shiny", package = "CrossICC"))
}

#' Title Summary the CrossICC-returned list to produce human-readable output
#'
#' @param result list-type CrossICC's return value.
#' @param iteration the result from which iteration time should be exhibited. Default is the last iteration.
#'
#' @return list contains: a matrix of genesets mapping to genes; a named atomic vetor of samples mapping to super clusters.
#' @export
#'
#' @examples
#' \donttest{
#' CrossICC.object <- CrossICC(example.matrices, max.iter = 20)
#' summary.CrossICC(CrossICC.object)
#' }
summary.CrossICC <- function(result) {

  temp.object<-result$clusters[[1]]
  # Get final gene2cluster set
  union.geneset<-result$unioned.genesets
  colnames(union.geneset)=c("Cluster","Genes")
  final.geneset<-union.geneset
  # Get final clusterSample result
  if(class(temp.object)=='list'){
    names(temp.object)=c()
    final.cluster<-do.call(c,temp.object)
  }else{
    final.cluster<-temp.object
  }

  data.matrx.list<-result$platforms

  colnames(final.geneset)=c("Cluster","Genes")
  list(gene.signatures = data.frame(final.geneset),
       clusters = final.cluster,normalized.matrix = data.matrx.list,order.gene=result$gene.order)
}
