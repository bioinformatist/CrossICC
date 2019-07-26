#' Run the default shiny report independently
#'
#' \code{runShinyCrossICC} run shiny independently
#' @author Qi Zhao
#' @export
#'
#' @example
#' \donttest{
#' runShinyCrossICC()
#' }
runShinyCrossICC <- function() {
  pkg.suggested <- c('rmarkdown', 'knitr', 'shiny', 'shinydashboard', 'shinyWidgets', "shinycssloaders", 'DT', 'ggthemes', 'ggplot2', 'pheatmap', 'RColorBrewer', 'tibble')
  checkPackages <- function(pkg){
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package pkg needed for shiny app. Please install it.",
           call. = FALSE)
    }
  }
  lapply(pkg.suggested, checkPackages)
  shiny::runApp(system.file("shiny", package = "CrossICC"))
}
