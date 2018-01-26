sil.width <- function(cc.matrix, k){
  # Make it silenced for finer strategy may produce non-square matrix
  dist.M <- suppressWarnings(as.dist(1 - cc.matrix))
  H <- hclust(dist.M, method = 'average')
  HC <- cutree(H, k = k)
  si <- cluster::silhouette(HC, dist.M)
  return(list(si, HC))
}
