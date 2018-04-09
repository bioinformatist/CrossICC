m.f.s <- function(platforms.list, filter.cutoff = 0.5, fdr.cutoff = 0.1, perform.mad = TRUE) {

  # Merge multiple probes for one gene here
  non.duplicates <- lapply(platforms.list, merge.duplicates)

  # no.same <- lapply(non.duplicates, function(x) x[apply(x[,-1], 1, function(y) !all(y==0)),])
  # no.same <- lapply(non.duplicates, function(x) x[apply(x[,-1], 1, function(y) !zero_range(x)),])

  # Remove rows with all values are same
  no.same <- lapply(non.duplicates, remove.all.same)

  genes.com <- com.feature(lapply(no.same, rownames), method = "overlap")

  if (length(platforms.list) > 1) {
    # Must pre-assigned here (deep copy, and must not be slice of list), or
    # mergeExprs will raise error (exactly its bug):
    # https://github.com/Bioconductor-mirror/MergeMaid/blob/baf0cfc0d370917d55c4b4adbc1d75c1141a3661/R/MergeMaid.R#L323
    # Error in dimnames(x) <- dn : length of 'dimnames' [2] not equal to array
    # extent MergeMaid.R has been fine-adjusted in this package.

    genes.com.list <- lapply(no.same, function(x) unlist(x)[genes.com, ])
    merged <- mergeExprs(genes.com.list)

    # fig.size <- (length(platforms.list) + 1) * 400
    # tiff('gene.cor.matrix.tiff', compression = 'lzw', res = 300, width = fig.size,
    # height = fig.size)

    graphics.off()
    pdf(NULL)
    par(mar=c(1,1,1,1))
    dev.control("enable")  # enable display list
    null.ic <- unlist(MergeMaid::intcorDens(merged))
    dev.off()

    pval.genes <- apply(as.matrix(rowMeans(MergeMaid::pairwise.cors(MergeMaid::intCor(merged,
      exact = FALSE)))), 1, function(x) {
      pval.cal(x, d = null.ic, alt = "g")
    })
    fdr.genes <- p.adjust(pval.genes, method = "fdr")
    genes.com.fdr <- names(which(fdr.genes < fdr.cutoff))

    if (length(genes.com.fdr) == 0) {
      stop("The dataset has no common signature gene which can pass FDR filtering!\n
           You may check if fake gene symbol/ID or wrongly annotated? Or try smaller FDR cutoff.")
    }

    if (perform.mad) {
      filter.genes <- lapply(no.same, function(x) filter.mad(x[genes.com.fdr,
                                                                ], p = filter.cutoff))
      filter.genes.com <- com.feature(filter.genes, method = "overlap")
    } else {
      filter.genes.com <- genes.com.fdr
    }

  } else {
    if (perform.mad) {
      filter.genes <- lapply(no.same, function(x) filter.mad(x, p = filter.cutoff))
      if(is.null(filter.genes)){
        warning("filter.cutoff is to large too get enough gene number in dataset")
      }
      filter.genes.com <- com.feature(filter.genes, method = "overlap")
    } else {
      filter.genes.com <- genes.com
    }
  }

  if (length(filter.genes.com) < 5) {
    stop("There are too few features in the dataset after filtering.\n
         Try to re-run with adjusted parameters or check your data source (too few overlapped features among matrices).")
  }

  filter.scale <- lapply(no.same, function(x) scale(t(scale(t(x[filter.genes.com,
    ])))))
  #filter outlier with range
  filter.scale<-lapply(filter.scale,function(x) replace(x, x < -2, -2) )
  filter.scale<-lapply(filter.scale,function(x) replace(x, x > 2, 2) )


  return(list(filtered.gene = filter.genes.com, filterd.scaled = filter.scale))
  # filter.genes.com
}
