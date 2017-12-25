balance.cluster <- function(xyz, cc, cluster.cutoff = 0.05, max.K = 6, plot = TRUE, iter){
  k <- vapply(cc, function(x) derive.clusternum(x, cluster.cutoff, max.K), 2333)

  xyz.k <- cor.cc(xyz, cc, k)

  silws <- unlist(lapply(2:max.K, function(x) mean(sil.width(xyz.k, x[[1]])[[1]][,3])))
  max.silw <- which.max(silws) + 1

  si <- sil.width(xyz.k, max.silw)
  hc <- si[[2]]

  hc.list <- list(x = hc[1:k[['x']]],
             y = hc[(k[['x']] + 1):(k[['x']]+k[['y']])],
             z = hc[(k[['x']] + k[['y']] + 1):(k[['x']] + k[['y']] + k[['z']])])

  cc.k.balanced <- cc.k.old <- lapply(names(k), function(x) cc[[x]][[k[[x]]]]$consensusClass)
  names(cc.k.balanced) <- names(cc.k.old) <- names(k)

  invisible(lapply(1:max.silw,
                   function(x) lapply(names(k),
                                      function(y) cc.k.balanced[[y]][which(cc.k.old[[y]] %in% which(hc.list[[y]] == x))] <- x)))
  if(plot){
    tiff(paste("cluster.centroid.correlation", iter, "tiff", sep = "."),
         width = 1600, height = 1600, res = 300, compression = 'lzw')
    gplots::heatmap.2(xyz.k, distfun = function(c) as.dist(1 - c),
              hclustfun = function(c) hclust(c, method = "average"),
              col = gplots::greenred, trace = "none", density.info = "none")
    dev.off()
    tiff(paste("silhouette.plot", iter, "tiff", sep = "."), res = 300, width = 1600, height = 1600)
    c1 <- rainbow(max.silw)
    plot(si[[1]], col = c1[hc])
    dev.off()
  }

  list(xyz.k, cc.k.balanced)
}
