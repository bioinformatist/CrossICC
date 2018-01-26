cor.cc <- function(xyz, cc, k, method = 'finer'){
  centroids.list <- lapply(names(k),
                           function(p) lapply(1:k[[p]],
                                              function(q) apply(xyz[[p]][, names(which(cc[[p]][[k[[p]]]]$consensusClass == q))],
                                                                1, median)))
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
