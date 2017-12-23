com.feature <- function(x, y, z, method = 'merge'){
  switch (method,
    'merge' = unique(c(x, y, z)),
    'overlap' = Reduce(intersect, list(x, y, z))
  )
}
