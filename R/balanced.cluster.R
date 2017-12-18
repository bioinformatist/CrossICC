balanced.cluster <- function(x, y, z, cc, cluster.cutoff = 0.05, max.K = 10, plot = TRUE, iter){
  k <- vapply(cc, function(x) derive.clusternum(x, cluster.cutoff, max.K), 2333)

}
