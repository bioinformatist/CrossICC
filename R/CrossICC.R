#' CrossICC: Automatically Aggregating and Summarizing Bioinformatics Results for Interactive Report.
#'
#'
#' @docType package
#' @name CrossICC
#' @import data.table MASS MergeMaid grDevices stats Biobase
#' @importFrom utils head
#' @importFrom methods new is
#' @importFrom graphics abline legend lines mtext pairs par text
#' @importFrom Biobase exprs
#' @importFrom dplyr count
NULL

#' The Main Function of the package
#'
#' @param output.dir the results' output directory.
#' @param max.iter the maximum number of iterations.
#' @param rep.runs number of subsamples during clustering.
#' @param pItem proportion of items to sample during clustering.
#' @param pFeature proportion of features to sample during clustering.
#' @param clusterAlg cluster algorithm. Could be 'hc' heirarchical (hclust),
#' 'pam' for paritioning around medoids, 'km' for k-means upon data matrix,
#' 'kmdist' for k-means upon distance matrices (former km option),
#' or a function that returns a clustering.
#' @param distance Could be 'pearson': (1 - Pearson correlation),
#' 'spearman' (1 - Spearman correlation),
#' 'euclidean', 'canberra', 'minkowski" or custom distance function.
#' @param cc.seed sets random seed for reproducible results.
#' @param filter.cutoff low variability (median absolute deviation (MAD)) cutoff threshold, default is 0.5.
#' @param fdr.cutoff cutoff value during fdr filtering.
#' @param cluster.cutoff cutoff value during determining cluster numbers.
#' @param ebayes.cutoff p-value cutoff when select differentially expressed probes.
#' @param study.names a vector containing all study names
#' @param ... all datasets (matrices is better)
#' @param ebayes.mode 'up' or 'both'. Choose only up-regulated genes or all differentially expressed genes when determining MDEGs. default is 'up'
#' @param use.shiny if TRUE, a shiny app will appear after running this main function. Note: You must keep output.dir with default value '~' for using shiny app.
#' @param cross object type when determining meta-cluster. Could be "cluster" for clusters by ConsencusClusterPlus, "sample" for samples or "none" (only used for single dataset).
#' @param max.K the maximum cluster number of ConsensusClusterPlus. Default is 10, but was set as number of samples when there're less than 10 samples.
#' @param skip.merge.dup skip merge multiple probes for one gene (duplicates) or not. Default is TRUE (it is highly recommended that user has their data pre-processed well).
#' @param skip.mm skip MergeMaid processing or not. Default is FALSE (not skip).
#' @param sil.filter silhouetee width filtering mode. Could be "soft" or "hard". If "hard", all negtive silhouetee width value will be set to 0. Default is "soft" (to do nothing).
#' @param heatmap.order gene order for heatmaps. Default is "up.based", with which genes will be arranged as up-regulated order in meta-clusters across all matrices. Or can be set to "concordant" for all in same order.
#' @param n.platform to filter the signature with it's meta-cluster group in platforms.
#' That is, if the parameter is set to 2 (default),
#' the signature (like hgnc symbol ESR1) in a certain meta-cluser (like K1) must exists more than 2 times among data of all platforms;
#' otherwise, it will not be reported.
#' @param skip.mfs by default, the datasets will be normalized at the start,
#' and the genes or features that have no or few contributions to the final clusters will be filtered out. To skip this process, you can set this parameter to TRUE. Only try when you're sure that you're working with pre-processed datasets.
#' @param com.mode mode for choose common features when pre-processing data. Could be "overlap" (use intersection, default) or "merge" (keep all features).
#' @param supercluster.method method for super-clustering. Default is 'hclust', can also be 'kmeans'.
#' @param overwrite if user allow overwrite result file? Default is FALSE.
#'
#' @return A nested list with iteration time as its name and list containing consensus cluster,
#' gene signature and balanced cluster as its value.
#'
#' @export
#'
#' @seealso \code{\link[ConsensusClusterPlus]{ConsensusClusterPlus}}
#'
#' @examples
#' data("demo.platforms")
#' CrossICC.obj <- CrossICC(demo.platforms, skip.mfs = TRUE, max.iter = 1, overwrite = TRUE)
CrossICC <- function(..., study.names, filter.cutoff = 0.5, fdr.cutoff = 0.001, output.dir = '~', max.K = 10, max.iter = 20, rep.runs = 1000, n.platform = 2,
                     pItem = 0.8, pFeature = 1, clusterAlg = "hc", distance = "euclidean", sil.filter = 'soft', heatmap.order = 'up.based', com.mode = 'overlap',
                     cc.seed = NULL, cluster.cutoff = 0.05, ebayes.cutoff = 0.1, ebayes.mode = 'up', cross = 'cluster', supercluster.method = 'hclust', skip.merge.dup = TRUE,
                     skip.mm = FALSE, skip.mfs = FALSE, use.shiny = FALSE, overwrite = FALSE){


  if (max.iter < 2) warning('Result from less than 2 times iteration may not make sense at all!')

  dir.create(output.dir, showWarnings = FALSE)


  # Check parameter values (sometimes users bring spelling mistake here)
  cross <- match.arg(cross, c("cluster", "sample", "none"))
  heatmap.order <- match.arg(heatmap.order, c('up.based', 'concordant'))
  # Get arguments as data.table
  arg.list <- unlist(as.list(match.call())[-1])

  graphics.off()

  arg <- list(...)
  result <- list()

  # It is not allowed just a matrix as input, for MergeMaid will raise an error:
  # Error in .dens.mergeExpressionSet(x = x, method = method, ...) :
  #   Number of studies in the mergeExpressionSet should not less than 2.
  # Already can auto-skip this filtering when there's only one sample.
  # if ((length(arg) == 1) & (is.element(class(arg[[1]]),"matrix"))) {
  #   stop("Number of studies should not less than 2.")
  if ((length(arg) == 1) & (is(arg[[1]], "list"))) {
    platforms.list <- unlist(arg, recursive = FALSE)
    if (length(platforms.list) == 1) {
      message('Only one matrix detected. MergeMaid will not work. Will skip cross analysis.')
      cross <- 'none'
    }
  } else {
    platforms.list <- lapply(arg, check.eSet)
  }

  # Merge, filter and scale here
  if (skip.mfs) {
    message('Merging, filtering and scaling are skipped.')
    platforms <- platforms.list
    filter.sig <- rownames(platforms[[1]])
  } else {
    cat(paste(date(), '--', 'Pre-processing data'), '\n')
    mfs.list <- m.f.s(platforms.list, fdr.cutoff = fdr.cutoff, filter.cutoff = filter.cutoff, skip.merge.dup = skip.merge.dup, skip.mm = skip.mm, com.mode = com.mode)
    platforms <- mfs.list$filterd.scaled
    filter.sig <- mfs.list$filter.sig
  }

  # If study.names is not defined or seems not OK, use automatically generated ones instead
  if (missing(study.names) || !is(study.names, "character") || length(study.names) != length(platforms)) {
    message('No study names provided or something goes wrong with your study names. Will use auto-generated study names instead.')
    study.names <- vapply(1:length(platforms), function(x) paste0('Matrix.', x), "Yu Fat is handsome")
  }

  names(platforms) <- study.names

  # Use the first datasets' feature names as default
  gene.sig <- rownames(platforms[[1]])

  # Use matrix save iteration time with signature number
  iter.sig <- data.table(Iteration = numeric(), Signatures = numeric())

  # To determine max k for ConsensusClusterPlus:
  # https://github.com/renzhonglu/ConsensusClusterPlus/blob/master/R/ConsensusClusterPlus.R#L514
  sampleN <- floor(min(vapply(platforms, function(x) dim(x)[2], 23333))*pItem)
  if (sampleN < 10) {
    max.K <- sampleN
  }

  iteration <- 1
  while(iteration <= max.iter){
    cat(paste(date(), iteration, sep=" -- start iteration: "), '\n')
    # Here a named vector is needed (for dir names), so (v)applys is necessary
    run.dir <- vapply(names(platforms),
                      function(x) paste(x, iteration, sep = "."),
                      "Yu Fat is handsome")

    cc <- vapply(names(platforms),
                 function(x) list(suppressMessages(ConsensusClusterPlus::ConsensusClusterPlus(platforms[[x]][gene.sig,],
                                                                                              maxK = max.K, reps=rep.runs, pItem=pItem,
                                                                                              pFeature=pFeature, title=run.dir[x],
                                                                                              clusterAlg=clusterAlg, distance=distance,
                                                                                              seed=cc.seed, plot = NULL))),
                 # ConsensusClusterPlus returns a list of 7 elements, but we need a nested list
                 list(rep(list('aString'), 7)))

    all.sig <- lapply(names(platforms), function(x) platforms[[x]][gene.sig,])
    names(all.sig) <- names(platforms)

    balanced.cluster <- balance.cluster(all.sig,
                                        cc = cc, cluster.cutoff = cluster.cutoff, sil.filter = sil.filter,
                                        max.K = max.K, cross = cross, supercluster.method = supercluster.method)

    ebayes.result <- lapply(names(all.sig),
                            function(x) ebayes(all.sig[[x]],
                                               switch(cross,
                                                      'cluster' = balanced.cluster[[1]][[x]],
                                                      'none' = unlist(balanced.cluster),
                                                      'sample' = balanced.cluster[[1]]),
                                               cutoff = ebayes.cutoff,
                                               mode = ebayes.mode))

    # Remove possible NULL element produced during ebayes step
    ebayes.result <- Filter(Negate(is.null), ebayes.result)

    gene.sig.all <- lapply(ebayes.result, function(x) rownames(x$full.m))
    geneset2gene <- lapply(ebayes.result, function(x) x$geneset2gene)
    # To avoid rbind raise "names do not match previous names"
    geneset2gene <- Filter(Negate(is.null), geneset2gene)
    geneset2gene <- lapply(geneset2gene , setNames , nm = c('super.cluster', 'signatures'))

    if (length(platforms) >= 2) {
      unioned.genesets <- data.table(do.call(rbind, geneset2gene))[, .(count = .N), by = 'super.cluster,signatures'][count >= n.platform,][, 1:2]
    } else {
      unioned.genesets <- as.matrix(unique(data.table(do.call(rbind, geneset2gene))))
    }

    pre.gene.sig <- gene.sig
    gene.sig <- com.feature(unlist(gene.sig.all), method = 'merge')

    merged.matrix <- do.call(cbind, all.sig)

    if (heatmap.order == 'concordant') {
      merged.cluster <- switch (cross,
                                'cluster' = do.call(c, balanced.cluster[[1]]),
                                'none' = unlist(balanced.cluster),
                                'sample' = balanced.cluster[[1]]
      )

      design <- model.matrix(~ 0 + factor(merged.cluster))
      k <- length(unique(merged.cluster))
      colnames(design) <- paste("K", seq_len(k), sep = "")
      K <- colnames(design)
      fit <- limma::lmFit(merged.matrix, design)

      x <- c()
      for(i in 1:(k - 1)){
        for(j in (i + 1):k){
          x <- c(x, paste(K[j], K[i], sep = "-"))
        }
      }

      for(i in seq_len(k)){
        x <- c(x, paste(K[i], paste("(", paste(K[-i], collapse="+"), ")", "/", k-1, sep=""), sep="-"))
      }

      contrast.matrix <- limma::makeContrasts(contrasts = x,levels = design)
      fit2 <- limma::contrasts.fit(fit,contrast.matrix)
      fit2 <- limma::eBayes(fit2)
      FC.matrix <- limma::topTable(fit2, number = 20000, adjust.method = 'BH')
      gene.order <- rownames(FC.matrix[order(-FC.matrix[,1]),])
      # gene.order <- sorted.genes[sorted.genes %in% gene.sig]
    } else {
      ebayes.result.2 <- lapply(names(all.sig),
                              function(x) ebayes(all.sig[[x]],
                                                 switch(cross,
                                                        'cluster' = balanced.cluster[[1]][[x]],
                                                        'none' = unlist(balanced.cluster),
                                                        'sample' = balanced.cluster[[1]]),
                                                 cutoff = 1,
                                                 mode = ebayes.mode))
      # Filtering: only keep 1 to all
      FC.tables <- lapply(ebayes.result.2, function(x) x$full.m[,grep('\\.\\.', names(x$full.m))])
      # Pre-ordering: From left to right, according to all columns
      sorted.tables <- lapply(FC.tables, function(x) x[do.call(order, lapply(1:NCOL(x), function(i) -x[, i])), ])
      gene.order <- lapply(sorted.tables, function(h) unique(unlist(lapply(1:ncol(h), function(x) row.names(h[h[,x] > 0,])))))
      # gene.order <- com.feature(unlist(lapply(sorted.genes, function(x) x[x %in% gene.sig])), method = 'merge')

    }


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
    er = ebayes.result,  # For test only
    arg.list = arg.list,  # named vector object of arguments
    platforms = platforms,  # For heatmap in shiny
    gene.signature = gene.sig,
    filter.sig = filter.sig,
    iter.sig = iter.sig,
    gene.order = gene.order,  # Sorted gene names by Fold-Change value, for heatmaps use
    clusters = balanced.cluster,
    geneset2gene = geneset2gene,
    unioned.genesets = as.matrix(unioned.genesets)
  )

  if (!overwrite & file.exists(file.path(output.dir, 'CrossICC.object.rds'))) {
    stop('Result file already existed!')
  }

  saveRDS(result, file = file.path(output.dir, 'CrossICC.object.rds'))
  cat("A CrossICC.object.rds file will be generated in home directory by default.
      Note that the previous file will be overridden.\n")

  cat(paste(date(), iteration - 1, sep=" -- Iteration finished! Iteration time for reaching convergence/limit: "), '\n')

  if (use.shiny) {
    # warning('Result list would not be returned when use.shiny = TRUE.')
    pkg.suggested <- c('ggalluvial', 'rmarkdown', 'knitr', 'shiny', 'shinydashboard', 'shinyWidgets', "shinycssloaders", 'DT', 'ggthemes', 'ggplot2', 'pheatmap', 'RColorBrewer', 'tibble')
    checkPackages <- function(pkg){
      if (!requireNamespace(pkg, quietly = TRUE)) {
        stop("Package ", pkg, " needed for shiny app. Please install it. All package needed are: ", paste(pkg.suggested, collapse = ' '), call. = FALSE)
      }
    }
    lapply(pkg.suggested, checkPackages)
    shiny::runApp(system.file("shiny", package = "CrossICC"))
  }
  result
}

#' Summary the CrossICC-returned list to produce human-readable output
#'
#' @param result list-type CrossICC's return value.
#'
#' @return list contains: a matrix of genesets mapping to genes; a named atomic vetor of samples mapping to super clusters.
#' @export
#'
#' @examples
#' data("demo.platforms")
#' CrossICC.object <- CrossICC(demo.platforms, skip.mfs = TRUE, max.iter = 1, overwrite = TRUE)
#' CrossICC.summary <- summaryCrossICC(CrossICC.object)
summaryCrossICC <- function(result) {

  temp.object<-result$clusters[[1]]
  # Get final gene2cluster set
  union.geneset<-result$unioned.genesets
  colnames(union.geneset)=c("Cluster","Genes")
  final.geneset<-union.geneset
  # Get final clusterSample result
  if(is(temp.object, 'list')){
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


#' Read file into CrossICC input
#'
#' @param files a list for filenames, usually a returned value of list.files() function
#' @return list contains matrices from each platform parsing from file provided .
#' @export
#'
#' @examples
#' files <- list.files(path=".", pattern = '.csv')
#' CrossICC.input <- CrossICCInput(files)
CrossICCInput <- function(files){
  if (is.character(files)) {
    dfinput <- function(x){
      outputdf <- as.matrix(data.frame(fread(x, stringsAsFactors = FALSE,  check.names = FALSE), row.names = 1, check.names = FALSE))
      return(outputdf)
    }
    testData <- lapply(files, dfinput)
    return(testData)
  }
  else {
    stop("Please ensure character type filenames were provided!")
  }
}

