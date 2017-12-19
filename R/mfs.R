mfs <- function(x, y, z, filter.cutoff = 0.5, fdr.cutoff = 0.001){
  invisible(sapply(list(x, y, z), function(x) setnames(x, 1, 'ID_REF')))
  genes.com <- com.feature(x[,1], y[,1], z[,1], method = 'overlap')
  merged <- MergeMaid::mergeExprs(as.matrix(subset(x, ID_REF %in% genes.com)[,-1]), as.matrix(subset(y, ID_REF %in% genes.com)[,-1]), as.matrix(subset(z, ID_REF %in% genes.com)[,-1]))
}
