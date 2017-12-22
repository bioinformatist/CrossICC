sil.width <- function(cc.matrix, k){
  dist.M <- as.dist(1 - cc.matrix)
  H <- hclust(dist.M, method = 'average')
  HC <- cutree(H, k = k)
  si <- cluster::silhouette(HC, dist.M)
  return(list(si, HC))
}
