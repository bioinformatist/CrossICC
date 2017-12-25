m.f.s <- function(x, y, z, filter.cutoff = 0.5, fdr.cutoff = 0.001){

  # Check if all input matrix is normal
  list.m <- lapply(list(x, y, z), check.eSet)

  # Merge multiple probes for one gene here
  non.duplicates <- lapply(list.m, merge.duplicates)

  # Add jitter for extremly same expressed features through samples (if existed)
  jittered <- lapply(non.duplicates, add.jitter)

  genes.com <- com.feature(rownames(jittered[[1]]),
                           rownames(jittered[[2]]),
                           rownames(jittered[[3]]), method="overlap")

  # Must pre-assigned here (deep copy, and must not be slice of list),
  # or mergeExprs will raise error (seems its bug):
  # Error in dimnames(x) <- dn :
  #   length of 'dimnames' [2] not equal to array extent
  x.genes.com <- jittered[[1]][genes.com,]
  y.genes.com <- jittered[[2]][genes.com,]
  z.genes.com <- jittered[[3]][genes.com,]

  merged <- MergeMaid::mergeExprs(x.genes.com, y.genes.com, z.genes.com)

  tiff("gene.cor.matrix.tiff", compression = 'lzw', res = 300, width = 1600, height = 1600)
  null.ic <- unlist(MergeMaid::intcorDens(merged))
  dev.off()

  pval.genes <- apply(as.matrix(rowMeans(MergeMaid::pairwise.cors(MergeMaid::intCor(merged, exact=FALSE)))), 1,
                      function(x){pval.cal(x,d=null.ic,alt="g")})
  fdr.genes <- p.adjust(pval.genes, method="fdr")
  genes.com.fdr <- names(which(fdr.genes < fdr.cutoff))
  filter.genes <- lapply(jittered, filter.mad, p = filter.cutoff)

  filter.genes.com <- com.feature(filter.genes[[1]],
                                  filter.genes[[2]],
                                  filter.genes[[3]], method="overlap")

  filter.scale <- lapply(jittered, function(x) scale(t(scale(t(x[filter.genes.com,])))))

  return(list(filtered.gene = filter.genes.com, filterd.scaled = filter.scale))
}
