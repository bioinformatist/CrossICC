balanced.cluster <- function(xyz, cc, cluster.cutoff = 0.05, max.K = 6, plot = TRUE, iter){
  k <- vapply(cc, function(x) derive.clusternum(x, cluster.cutoff, max.K), 2333)

  xyz.k <- cor.cc(xyz, cc, k)

  silws <- unlist(lapply(2:max.K, function(x) mean(sil.width(xyz.k, x[[1]])[[1]][,3])))
  max.silw <- which.max(silws) + 1

  si <- sil.width(xyz.k, max.silw)
  hc <- si[[2]]

  hc <- list(x = hc[1:k[['x']]],
             y = hc[(k[['x']]+1):(k[['x']]+k[['y']])],
             z = hc[(k[['x']]+k[['y']]+1):(k[['x']]+k[['y']]+k[['z']])])

  cc.k.balanced <- cc.k.old <- lapply(names(k), function(x) cc[[x]][[k[[x]]]]$consensusClass)
  names(cc.k.balanced) <- names(cc.k.old) <- names(k)

  invisible(lapply(1:max.silw,
                   function(x) lapply(names(k),
                                      function(y) cc.k.balanced[[y]][which(cc.k.old[[y]] %in% which(hc[[y]] == x))] <- x)))
}
