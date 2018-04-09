balance.cluster <- function(sig.list, cc, cluster.cutoff = 0.05, max.K = NULL, plot = TRUE, iter, method){
  k <- vapply(cc, function(x) derive.clusternum(x, cluster.cutoff), 2333)

  # Max cluster number must be refined here, for silhouette statistics are only defined if 2 <= k <= n-1.
  # Here, n is sum(k), Qi Fat suggests that max.K should be as half as n by default
  if (is.null(max.K)) {
    max.K <- ceiling(sum(k) / 2)
  } else {
    if ((sum(k) - 1) < max.K) {
      warning("User defined max super cluster number is larger than real value, will refine it to (n - 1)!")
      max.K <- sum(k) - 1
    }
  }

  if (max.K < 2) {
    max.K <- 2
  }

  all.k <- cor.cc(sig.list, cc, k, method = method)

  silws <- unlist(lapply(2:max.K, function(x) mean(sil.width(all.k, x)[[1]][,3])))
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

    # win.metafile()
    pdf(NULL)
    dev.control('enable') # enable display list
    c1 <- rainbow(max.silw)
    plot(si[[1]], col = c1[hc])
    silhouette <- recordPlot()
    dev.off()
    # tiff(paste("silhouette.plot", iter, "tiff", sep = "."), res = 300, width = 1600, height = 1600)
  }

  list(  # all.k,
       clusters = cc.k.balanced, heatmap = heatmap, silhouette = silhouette)
}

derive.clusternum <- function(consencus.result, cutoff = 0.05, maxK = 7){
  for(k in 2:maxK){
    A <- cal.auc(consencus.result, k)
    if(k == 2){
      delta.A <- A
      pre.A <- A
    } else delta.A <- (A-pre.A) / pre.A
    if (delta.A <= cutoff) {
      # Return value must be larger than 2!
      if ((k - 1) < 2) {
        return(2)
      } else {
        return(k - 1)
      }
    } else if (k == maxK){
      return(maxK)
    }
    pre.A <- A
  }
}

cal.auc <- function(consencus.result, k){
  consencus.matrix <- consencus.result[[k]]$consensusMatrix
  cdf <- ecdf(c(consencus.matrix))

  # Get diagonal after removing first column and last row
  t <- sort(diag(consencus.matrix[-nrow(consencus.matrix), -1]))

  # Make a "lagged" t for performing vectorized calculation
  lagged.t <- c(NA, head(t, -1))
  # cdf object is closure that is not subsettable, should be called once used
  cdf.t <- cdf(t)
  sum((t[-1] - lagged.t[-1]) * cdf.t[-1])
}

cor.cc <- function(xyz, cc, k, method = 'finer'){
  # cat(xyz$Matrix.1, "\n")
  # cat(names(k), k, "\n")
  # cat(xyz$Matrix.1[, names(which(cc$Matrix.1[[5]]$consensusClass == 1))], "\n")
  centroids.list <- lapply(names(k),
                           function(p) lapply(1:k[[p]],
                                              # Thanks to https://stackoverflow.com/questions/28423275/dimx-must-have-a-positive-length-when-applying-function-in-data-frame/28423503#28423503
                                              function(q) apply(xyz[[p]][, names(which(cc[[p]][[k[[p]]]]$consensusClass == q)), drop = FALSE],
                                                                1, median, na.rm = TRUE)))
  centroids <- as.data.table(unlist(centroids.list,
                                    recursive = FALSE))

  centroids.names.list <- lapply(names(k),
                                 function(x) lapply(1:k[[x]],
                                                    function(y) paste(x, k[[x]], y, sep = '.')))
  centroids.names <- unlist(centroids.names.list)
  setnames(centroids, centroids.names)

  switch (method,
          'balanced' = cor(centroids, method="pearson"),
          'finer' = cor(centroids, do.call(cbind, xyz))
  )
}

sil.width <- function(cc.matrix, k){
  # Make it silenced for finer strategy may produce non-square matrix
  dist.M <- suppressWarnings(as.dist(1 - cc.matrix))
  H <- hclust(dist.M, method = 'average')
  HC <- cutree(H, k = k)
  si <- cluster::silhouette(HC, dist.M)
  return(list(si, HC))

  dist.M <- suppressWarnings(dist(cc.matrix))
  H <- hclust(dist.M, method = 'average')
  HC <- cutree(H, k = k)
  si <- cluster::silhouette(HC, dist.M)
  return(list(si, HC))
}
