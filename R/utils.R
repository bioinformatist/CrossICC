check.eSet <- function(x){
  if(!is.numeric(x[[1]])){
    stop("The first element is not numeric:\n
         Currently, we only support matrix (data.frame) with sample names as column names
         and feature names as row names.")
  }
  if(!is.matrix(x)){
    as.matrix(x)
  } else {
    x
  }
}

get.max.var.row <- function(x){
  # Thanks to https://stackoverflow.com/users/3001626/david-arenburg
  # at https://stackoverflow.com/a/25100036 for this line.
  rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
  # TODO: return row with maximum var, then embed it into below
}

merge.duplicates <- function(x, method = "median") {
  switch(method,
         "median" = aggregate(x, by=list(row.names(x)), FUN=median),
         "mean" = aggregate(x, by=list(row.names(x)), FUN=median),
         "max" = aggregate(x, by=list(row.names(x)), FUN=max),
         "min" = aggregate(x, by=list(row.names(x)), FUN=min),
         "variance" = )
}
