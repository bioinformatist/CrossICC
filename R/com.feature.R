#' Title Extract Common Features
#'
#' @param x a vector, i.e., should be x[[1]] if x is a column of data.table/data.frame.
#' @param y ditto.
#' @param z ditto.
#' @param method should be "merge" or "overlap".
#'
#' @return
com.feature <- function(x, y, z, method = 'merge'){
  switch (method,
    'merge' = unique(c(x, y, z)),
    'overlap' = Reduce(intersect, list(x, y, z))
  )
}
