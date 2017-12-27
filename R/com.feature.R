# ... stand for different atomic vectors
com.feature <- function(..., method = 'merge'){
  switch (method,
    'merge' = unique(c(...)),
    'overlap' = Reduce(intersect, ...)
  )
}
