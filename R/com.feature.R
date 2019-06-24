# ... stand for different atomic vectors
com.feature <- function(..., method = 'merge'){
  switch(method,
    'merge' = unique(c(...)),  # input should be a vector.
    'overlap' = Reduce(intersect, ...)
  )
}
