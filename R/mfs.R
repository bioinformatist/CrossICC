mfs <- function(x, y, z, filter.cutoff = 0.5, fdr.cutoff = 0.001){

  # Check if all input matrix is normal
  invisible(lapply(list(x, y, z), check.eSet))

  # TODO: Merge multiple probes for one gene here

  genes.com <- com.feature(rownames(x), rownames(y), rownames(z), method="overlap")

  # Must pre-assigned here (deep copy), or mergeExprs will raise error (seems its bug):
  # Error in dimnames(x) <- dn :
  #   length of 'dimnames' [2] not equal to array extent
  x.genes.com <- x[genes.com,]
  y.genes.com <- y[genes.com,]
  z.genes.com <- z[genes.com,]

  merged <- mergeExprs(x.genes.com,y.genes.com,z.genes.com)

}
