# x is a list object only containing one element (used in apply)
check.eSet <- function(x) {
  if (!is.numeric(x[1])) {
    stop("The first element from your data is not numeric:\n
         Currently, we only support matrix (data.frame) with sample names as column names
         and feature names as row names.")
  } else {
    if (!is.matrix(x)) {
      as.matrix(x)
    } else {
      x
    }
  }
}

matrix.please <- function(x) {
  # Thanks to https://stackoverflow.com/a/25695207
  m<-as.matrix(x[,-1])
  rownames(m)<-x[,1]
  m
}

get.max.var.rows <- function(x) {
  DT <- data.table(row.names(x), x, var = apply(x, 1, var))
  # Thanks to https://stackoverflow.com/a/29497254
  DT <- setDT(DT)[, .SD[which.max(var)], by=V1]
  m <- matrix.please(DT)
  m[,-dim(m)[2]]
}

# Thanks to https://stackoverflow.com/a/4752580 for this function
zero_range <- function(x, tol = .Machine$double.eps ^ 0.5) {
  if (length(x) == 1) return(TRUE)
  x <- range(x) / mean(x)
  isTRUE(all.equal(x[1], x[2], tolerance = tol))
}

merge.duplicates <- function(x, method = "median") {
  switch(method,
         # Thanks to https://support.bioconductor.org/u/3556/
         # at https://support.bioconductor.org/p/41288/#41289 for this statement.
         "median" = matrix.please(aggregate(x, by=list(row.names(x)), FUN = median)),
         "mean" = matrix.please(aggregate(x, by=list(row.names(x)), FUN = mean)),
         "max" = matrix.please(aggregate(x, by=list(row.names(x)), FUN = max)),
         "min" = matrix.please(aggregate(x, by=list(row.names(x)), FUN = min)),
         "var" = get.max.var.rows(x))
}

add.jitter <- function(x) {
  t(apply(x, 1, function(x) if(zero_range(x) == TRUE){jitter(x)}else{x}))
}


remove.all.same <- function(m) {
  keep <- apply(m, 1, function(x) length(unique(x[!is.na(x)])) != 1)
  m[keep, ]
}
