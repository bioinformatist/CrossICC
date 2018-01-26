cluster.centroid <- function(x, gene.signature, cluster){
  # gene.signature <- intersect(rownames(x),gene.signature)
  x.scale <- t(scale(t(x[gene.signature,])))
  centroids <- c()
  for(i in unique(cluster)){
    centroids <- cbind(centroids,
                       apply(x.scale[gene.signature, names(which(cluster==i))], 1, median, na.rm=FALSE)
                       )
  }
  colnames(centroids) <- unique(cluster)
  return(centroids)
}
