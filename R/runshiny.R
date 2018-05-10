#' Run the default shiny report
#'
#' \code{runShinyCrossICC} run shiny independently
#' @author Qi Zhao
#' @export
runShinyCrossICC<-function(){
  shiny::runApp(system.file("shiny", package = "CrossICC"))
}
