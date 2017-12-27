m.f.s <- function(platforms.list, filter.cutoff = 0.5, fdr.cutoff = 0.001){

  # Merge multiple probes for one gene here
  non.duplicates <- lapply(platforms.list, merge.duplicates)

  # Add jitter for extremly same expressed features through samples (if existed)
  jittered <- lapply(non.duplicates, add.jitter)

  genes.com <- com.feature(lapply(jittered, rownames), method = 'overlap')

  genes.com.list <- lapply(jittered, function(x) unlist(x)[genes.com,])

  # Must pre-assigned here (deep copy, and must not be slice of list),
  # or mergeExprs will raise error (exactly its bug):
  # https://github.com/Bioconductor-mirror/MergeMaid/blob/baf0cfc0d370917d55c4b4adbc1d75c1141a3661/R/MergeMaid.R#L323
  # Error in dimnames(x) <- dn :
  #   length of 'dimnames' [2] not equal to array extent
  # MergeMaid.R has been fine-adjusted in this package.

  merged <- MergeMaid::mergeExprs(x.genes.com, y.genes.com, z.genes.com)

  tiff("gene.cor.matrix.tiff", compression = 'lzw', res = 300, width = 1600, height = 1600)
  null.ic <- unlist(MergeMaid::intcorDens(merged))
  dev.off()

  pval.genes <- apply(as.matrix(rowMeans(MergeMaid::pairwise.cors(MergeMaid::intCor(merged, exact=FALSE)))), 1,
                      function(x){pval.cal(x, d=null.ic, alt="g")})
  fdr.genes <- p.adjust(pval.genes, method="fdr")
  genes.com.fdr <- names(which(fdr.genes < fdr.cutoff))
  filter.genes <- lapply(jittered, filter.mad, p = filter.cutoff)

  filter.genes.com <- com.feature(filter.genes[[1]],
                                  filter.genes[[2]],
                                  filter.genes[[3]], method="overlap")

  filter.scale <- lapply(jittered, function(x) scale(t(scale(t(x[filter.genes.com,])))))

  return(list(filtered.gene = filter.genes.com, filterd.scaled = filter.scale))
}
