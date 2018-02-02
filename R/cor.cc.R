cor.cc <- function(xyz, cc, k, method = 'finer'){
  # cat(xyz$Matrix.1, "\n")
  # cat(names(k), k, "\n")
  # cat(xyz$Matrix.1[, names(which(cc$Matrix.1[[5]]$consensusClass == 1))], "\n")
  centroids.list <- lapply(names(k),
                           function(p) lapply(1:k[[p]],
                                              # Thanks to https://stackoverflow.com/questions/28423275/dimx-must-have-a-positive-length-when-applying-function-in-data-frame/28423503#28423503
                                              function(q) apply(xyz[[p]][, names(which(cc[[p]][[k[[p]]]]$consensusClass == q)), drop = FALSE],
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
