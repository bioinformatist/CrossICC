filter.mad <- function(x, p = 0.5, method = 'absolute'){
  x.mad <- apply(x, 1, function(xr) mad(xr[!is.na(xr)]))
  x.mad.rank <- rank(-x.mad)
  switch (method,
    "percent" = names(x.mad.rank[x.mad.rank < p * length(x.mad)]),
    "absolute" = names(x.mad[x.mad > p])
  )
}
