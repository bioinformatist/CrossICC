options(shiny.maxRequestSize=1024*1024^2)
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# add annotation
#

suppressMessages(library(shiny))
suppressMessages(library(pheatmap))
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
  #animate bar ---
  output$interationNumberForplot <- renderUI({
    if(is.null(InterationResult())){
      tagList(div())
    }else{
      df=InterationResult()
      iter.num<-length(df)
      tagList(
        sliderInput("iterslided","Total Iteration Time",min=1,max=iter.num,value = 1,step=2,animate=T)
      )
    }
  })
  #heatmap control option of platform selection ui----
  output$expressionHeatmapSelectPlatform <- renderUI({
    if(is.na(InterationResult())){
      tagList(div())
    }else{
      df=InterationResult()
      platformnamelist<-names(df[[input$iterslided]]$heatmaps)
      tagList(
        selectInput("SelectPL", "SelectPlatform", choices=platformnamelist, selected = platformnamelist[1], multiple = FALSE)
      )
    }
  })

  #input data ----
  inputdata<-  reactive({
    example<- readRDS(file = path.expand('~/CrossICC.object.rds'))
    inFile <- input$file1
    CrossICC.object<-NULL
    if (!is.null(inFile)){
      CrossICC.object<-readRDS(file = inFile$datapath)
    }
    switch(input$dataset,
           "default" = example,
           "upload" = CrossICC.object
    )
  })
  #interation CrossICC
  InterationResult <- reactive({

    # Create a Progress object
    if(input$submit!=0)
      CrossICC.object=inputdata()


  })

  #summary crossICC result
  output$OutputResultSignature <- renderDataTable({
    validate(
      need(!is.null(InterationResult()), "Please upload a correct CrossICC output file in RDA format, which can be found at default output path of CrossICC function or user defined path.")
    )
    tempr<-InterationResult()
    len<-length(tempr)
    s<-summary.CrossICC(tempr)
    s$gene.signatures
    })
  output$OutputClusterResult <- renderPrint({
    validate(
      need(!is.null(InterationResult()), "Please upload a correct CrossICC output file in RDA format, which can be found at default output path of CrossICC function or user defined path.")
    )
    tempr<-InterationResult()
    len<-length(tempr)
    s<-summary.CrossICC(tempr)
    s$clusters
  })

  #Plot functions
  output$superclusterPlot<-renderPlot({
    validate(
      need(!is.null(InterationResult()), "Please upload a correct CrossICC output file in RDA format, which can be found at default output path of CrossICC function or user defined path.")
    )
    fuck<-InterationResult()
    grid::grid.newpage()
    grid::grid.draw(fuck[[input$iterslided]]$clusters$heatmap$gtable)

  })
  output$Silhouette<-renderPlot({
    validate(
      need(!is.null(InterationResult()), "Please upload a correct CrossICC output file in RDA format, which can be found at default output path of CrossICC function or user defined path.")
    )
    fuck<-InterationResult()

    fuck[[1]]$clusters$silhouette
  })
  output$clusterexpress<-renderPlot({
    validate(
      need(!is.null(InterationResult()), "Please upload a correct CrossICC output file in RDA format, which can be found at default output path of CrossICC function or user defined path.")
    )
    fuck<-InterationResult()
    grid::grid.newpage()
    grid::grid.draw(fuck[[input$iterslided]]$heatmaps[[input$SelectPL]]$gtable)
  })
  #Download functions
  output$DownloadSuperclusterPlot<-downloadHandler(
    filename = function() {
      paste("SuperclusterPlot_", Sys.time(), '.pdf', sep='')
    },

    content = function(file) {
      pdf(file)
      fuck<-InterationResult()
      grid::grid.newpage()
      grid::grid.draw(fuck[[input$iterslided]]$clusters$heatmap$gtable)
      dev.off()
    },
    contentType = 'image/pdf'
  )
  output$DownloadSilhouette<-downloadHandler(
    filename = function() {
      paste("Silhouette_", Sys.time(), '.pdf', sep='')
    },

    content = function(file) {
      pdf(file)
      fuck<-InterationResult()
     fuck[[input$iterslided]]$clusters$silhouette
      dev.off()
    },
    contentType = 'image/pdf'
  )
  output$DownloadClusterexpressPlot<-downloadHandler(
    filename = function() {
      paste("Clusterexpress_", Sys.time(), '.pdf', sep='')
    },

    content = function(file) {
      pdf(file)
      fuck<-InterationResult()
      grid::grid.newpage()
      grid::grid.draw(fuck[[input$iterslided]]$heatmaps[[input$SelectPL]]$gtable)
      dev.off()
    },
    contentType = 'image/pdf'
  )

})
