balance.cluster <- function(sig.list, cc, cluster.cutoff = 0.05, max.K = 6, plot = TRUE, iter){
  k <- vapply(cc, function(x) derive.clusternum(x, cluster.cutoff, max.K), 2333)

  if ((sum(k) - 1) < max.K) {
    max.K = sum(k) - 1
  }

  all.k <- cor.cc(sig.list, cc, k)

  silws <- unlist(lapply(2:max.K, function(x) mean(sil.width(all.k, x[[1]])[[1]][,3])))
  max.silw <- which.max(silws) + 1

  si <- sil.width(all.k, max.silw)
  hc <- si[[2]]

  hc.list <- lapply(names(k), function(x) hc[grep(x, names(hc))])
  names(hc.list) <- names(k)

  cc.k.balanced <- cc.k.old <- lapply(names(k), function(x) cc[[x]][[k[[x]]]]$consensusClass)
  names(cc.k.balanced) <- names(cc.k.old) <- names(k)

  invisible(lapply(1:max.silw,
                   function(x) lapply(names(k),
                                      # Must use <<- here for scope restrinction
                                      function(y) cc.k.balanced[[y]][which(cc.k.old[[y]] %in% which(hc.list[[y]] == x))] <<- x)))
  if(plot){
    # tiff(paste("cluster.centroid.correlation", iter, "tiff", sep = "."),
    #      width = 1600, height = 1600, res = 300, compression = 'lzw')
    # op <- par(mar = c(4, 2, 2, 3))
    # win.metafile()
    # dev.control('enable') # enable display list
    # gplots::heatmap.2(all.k, distfun = function(c) as.dist(1 - c),
    #                   hclustfun = function(c) hclust(c, method = "average"),
    #                   col = gplots::greenred, trace = "none", density.info = "none", margins = c(9, 9))
    heatmap <- pheatmap::pheatmap(all.k,
                                  border_color = NA,
                                  colorRampPalette(c("green", "black", "red"))(50))
    # heatmap <- recordPlot()
    # dev.off()
    # par(op)
    # replayPlot(obj)

    win.metafile()
    dev.control('enable') # enable display list
    c1 <- rainbow(max.silw)
    plot(si[[1]], col = c1[hc])
    silhouette <- recordPlot()
    dev.off()
    # tiff(paste("silhouette.plot", iter, "tiff", sep = "."), res = 300, width = 1600, height = 1600)
  }

  list(all.k, cc.k.balanced, heatmap = heatmap, silhouette = silhouette)
}
