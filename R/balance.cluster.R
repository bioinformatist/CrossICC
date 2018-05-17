balance.cluster <- function(sig.list, cc, cluster.cutoff = 0.05, max.K = NULL, cross = 'sample', sil.filter = 'soft'){
  k <- vapply(cc, function(x) derive.clusternum(x, cluster.cutoff, maxK = max.K), 2333)

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

  if (cross == 'none') {
    clusters <- lapply(names(k), function(x) cc[[x]][[k[[x]]]]$consensusClass)
    return(clusters)
  }

  all.k <- cor.cc(sig.list, cc, k, cross = cross)

  # Sometimes correlation matrix row number is less than max cutree k
  if (dim(all.k)[1] < max.K) max.K <- dim(all.k)[1]

  silws <- unlist(lapply(2:max.K, function(x) mean(sil.width(all.k, x, cross = cross)[[1]][,3])))
  max.silw <- which.max(silws) + 1

  si <- sil.width(all.k, max.silw, cross = cross)
  hc <- si[[2]]



  if (cross == "cluster") {
    hc.list <- lapply(names(k), function(x) hc[grep(x, names(hc))])
    names(hc.list) <- names(k)

    cc.k.balanced <- cc.k.old <- lapply(names(k), function(x) cc[[x]][[k[[x]]]]$consensusClass)
    names(cc.k.balanced) <- names(cc.k.old) <- names(k)

    invisible(lapply(1:max.silw,
                     function(x) lapply(names(k),
                                        # Must use <<- here for scope restrinction
                                        function(y) cc.k.balanced[[y]][which(cc.k.old[[y]] %in% which(hc.list[[y]] == x))] <<- x)))
  }
  list(clusters = switch(cross,
                         "cluster" = cc.k.balanced,
                         'sample' = hc
  ),
  all.k = all.k, # For one heatmap
  silhouette = si[[1]])
}

derive.clusternum <- function(consencus.result, cutoff = 0.05, maxK = 10){
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

cor.cc <- function(xyz, cc, k, cross = 'sample'){
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

  switch (cross,
          'cluster' = cor(centroids, method = "pearson"),
          'sample' = t(cor(centroids, do.call(cbind, xyz)))
  )
}

sil.width <- function(cc.matrix, k, cross = 'sample', sil.filter){
  dist.M <- switch (cross,
                    'cluster' = as.dist(1 - cc.matrix),
                    'sample' = dist(cc.matrix)
  )

  H <- hclust(dist.M, method = 'average')
  HC <- cutree(H, k = k)
  si <- cluster::silhouette(HC, dist.M)
  if (sil.filter == 'hard') {
    si[,3][si[,3] < 0.5] <- 0
  }
  return(list(si, HC))
}
