#' Title
#'
#' @param x
#' @param d
#' @param alt
#'
#' @return
pval.cal <- function(x, d, alt="g"){
  switch (alt,
    "g" = 1 - pnorm(x, mean(d), sd(d)),
    "l" = pnorm(x,mean(d),sd(d))
  )
}
