
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# add annotation
#

suppressMessages(library(shiny))
suppressMessages(library(pheatmap))
# suppressMessages(library(CrossICC))
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
 #animate bar ---
  output$interationNumberForplot <- renderUI({
    df=InterationResult()
    iter.num<-length(df)
    tagList(
      sliderInput("iterslided","Total Iteration Time",min=1,max=iter.num,value = 1,step=2,animate=T)
    )
  })
#heatmap control option of platform selection ui----
  output$expressionHeatmapSelectPlatform <- renderUI({
    df=InterationResult()
    platformnamelist<-names(df[[input$iterslided]]$heatmaps)
    tagList(
      selectInput("SelectPL", "SelectPlatform", choices=platformnamelist, selected = platformnamelist[1], multiple = FALSE)
    )
  })

#input data ----
  inputdata<-  reactive({
    example.matrices
  })
  #interation CrossICC
  InterationResult <- reactive({
    if (input$submit != 0) {
      # Create a Progress object
      progress <- shiny::Progress$new()
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())

      progress$set(message = "Run iteration", value = 0)
      CrossICC.object <- CrossICCshiny(example.matrices, max.iter = 20,progress=progress)

      CrossICC.object
    } else{
      NULL
    }

  })

#Plot functions
  output$superclusterPlot<-renderPlot({
    fuck<-InterationResult()
    grid.newpage()
    grid.draw(fuck[[input$iterslided]]$balanced.cluster$heatmap$gtable)
  })
  output$Silhouette<-renderPlot({
    fuck<-InterationResult()
    replayPlot(fuck[[input$iterslided]]$balanced.cluster$silhouette)
  })
  output$clusterexpress<-renderPlot({
    fuck<-InterationResult()
    grid::grid.draw(fuck[[input$iterslided]]$heatmaps[[input$SelectPL]]$gtable)
  })


})
