
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# add annotation
#

suppressMessages(library(shiny))
suppressMessages(library(shinyjs))
suppressMessages(library(heatmaply))
suppressMessages(library(CrossICC))

shinyServer(function(session,input, output) {
# ui setting----
  output$workflowImage <- renderImage({
    return(list(
      src = "www/images/workflow.png",
      contentType = "image/png",
      height = 600,

      alt = "Face"
    ))
  })


#input data ----
inputdata<- reactive(
  example.matrices
)


})
