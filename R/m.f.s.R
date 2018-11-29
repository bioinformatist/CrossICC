m.f.s <- function(platforms.list, filter.cutoff = 0.5, fdr.cutoff = 0.1, perform.mad = TRUE, skip.merge.dup = FALSE, skip.mm = FALSE, com.mode = 'overlap') {

  set.seed(920304)

  # Check if has NAs in matrices
  if (!all(sapply(platforms.list, function(x) !any(is.na(x))))) {
    stop("Your list has as least one matrix contains NA(s)!\n
         You should process missing values first, for details, see https://github.com/bioinformatist/CrossICC#faqs")
  }

  if (!skip.merge.dup) {
    # Merge multiple probes for one gene here
    cat(paste(date(), '--', 'Merging multiple probes for one feature'), '\n')
    non.duplicates <- lapply(platforms.list, merge.duplicates)
  } else {
    non.duplicates <- platforms.list
  }

  # no.same <- lapply(non.duplicates, function(x) x[apply(x[,-1], 1, function(y) !all(y==0)),])
  # no.same <- lapply(non.duplicates, function(x) x[apply(x[,-1], 1, function(y) !zero_range(x)),])

  # Remove rows with all values are same
  cat(paste(date(), '--', 'Removing features with no variance'), '\n')
  no.same <- lapply(non.duplicates, remove.all.same)

  filter.sig <- rbind(vapply(platforms.list, function(x) dim(x)[1], 2333), vapply(non.duplicates, function(x) dim(x)[1], 2333), vapply(no.same, function(x) dim(x)[1], 2333))
  rownames(filter.sig) <- c('Original', 'Duplicates removed', 'No variance removed')

  genes.com <- com.feature(lapply(no.same, rownames), method = "overlap")

  filter.sig <- rbind(filter.sig, rep(length(genes.com), length(platforms.list)))
  rownames(filter.sig)[nrow(filter.sig)] <- 'Common signature'

  if ((length(platforms.list) > 1) & (skip.mm == FALSE)) {
    # Must pre-assigned here (deep copy, and must not be slice of list), or
    # mergeExprs will raise error (exactly its bug):
    # https://github.com/Bioconductor-mirror/MergeMaid/blob/baf0cfc0d370917d55c4b4adbc1d75c1141a3661/R/MergeMaid.R#L323
    # Error in dimnames(x) <- dn : length of 'dimnames' [2] not equal to array
    # extent MergeMaid.R has been fine-adjusted in this package.

    genes.com.list <- lapply(no.same, function(x) unlist(x)[genes.com, ])

    cat(paste(date(), '--', 'Performing MergeMaid'), '\n')
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

    filter.sig <- rbind(filter.sig, rep(length(genes.com.fdr), length(platforms.list)))
    rownames(filter.sig)[nrow(filter.sig)] <- 'After MergeMaid FDR filtering'

    if (length(genes.com.fdr) == 0) {
      stop("The dataset has no common signature gene which can pass FDR filtering!\n
           You may check if fake gene symbol/ID or wrongly annotated? Or try smaller FDR cutoff.")
    }

    if (perform.mad) {
      cat(paste(date(), '--', 'Performing MAD filtering'), '\n')
      filter.genes <- lapply(no.same, function(x) filter.mad(x[genes.com.fdr,
                                                               ], p = filter.cutoff))
      if (com.mode == 'overlap') {
        filter.genes.com <- com.feature(filter.genes, method = "overlap")
      } else {
        filter.genes.com <- com.feature(unlist(filter.genes), method = "merge")  # remove unlist function when using `overlap` parameter
      }

      filter.sig <- rbind(filter.sig, rep(length(filter.genes.com), length(platforms.list)))
      rownames(filter.sig)[nrow(filter.sig)] <- 'After MAD filtering'

    } else {
      filter.genes.com <- genes.com.fdr
    }

  } else {
    if (perform.mad) {
      cat(paste(date(), '--', 'Performing MAD filtering'), '\n')
      filter.genes <- lapply(no.same, function(x) filter.mad(x, p = filter.cutoff))
      if(is.null(filter.genes)){
        warning("filter.cutoff is to large too get enough gene number in dataset")
      }

      if (com.mode == 'overlap') {
        filter.genes.com <- com.feature(filter.genes, method = "overlap")
      } else {
        filter.genes.com <- com.feature(unlist(filter.genes), method = "merge")  # remove unlist function when using `overlap` parameter
      }

      filter.sig <- rbind(filter.sig, rep(length(filter.genes.com), length(platforms.list)))
      rownames(filter.sig)[nrow(filter.sig)] <- 'After MAD filtering'
    } else {
      filter.genes.com <- genes.com
    }
  }

  if (length(filter.genes.com) < 5) {
    stop("There are too few features in the dataset after filtering.\n
         Try to re-run with adjusted parameters or check your data source (too few overlapped features among matrices).")
  }

  cat(paste(date(), '--', 'Scaling'), '\n')
  filter.scale <- lapply(no.same, function(x) scale(t(scale(t(x[filter.genes.com,
    ])))))

  # Filter outlier with range
  filter.scale<-lapply(filter.scale,function(x) replace(x, x < -2, -2))
  filter.scale<-lapply(filter.scale,function(x) replace(x, x > 2, 2))

  filter.sig <- rbind(filter.sig, filter.sig[nrow(filter.sig),])
  rownames(filter.sig)[nrow(filter.sig)] <- 'After pre-processing / Before Iteration'

  list(filtered.gene = filter.genes.com, filterd.scaled = filter.scale, filter.sig = filter.sig)
  # filter.genes.com
}

pval.cal <- function(x, d, alt="g"){
  switch (alt,
          "g" = 1 - pnorm(x, mean(d), sd(d)),
          "l" = pnorm(x,mean(d),sd(d))
  )
}

filter.mad <- function(x, p = 0.5, method = 'absolute'){
  x.mad <- apply(x, 1, function(xr) mad(xr[!is.na(xr)]))
  x.mad.rank <- rank(-x.mad)
  switch (method,
          "percent" = names(x.mad.rank[x.mad.rank < p * length(x.mad)]),
          "absolute" = names(x.mad[x.mad > p])
  )
}
