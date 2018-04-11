options(shiny.maxRequestSize=1024*1024^2)
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# add annotation
#

suppressMessages(library(shiny))
suppressMessages(library(pheatmap))
suppressMessages(library(CrossICC))
suppressMessages(library(RColorBrewer))
shinyServer(function(session,input, output) {
# CrossICC panel functions ----
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
        sliderInput("iterslided","Total Iteration Time",min=1,max=iter.num,value = iter.num,step=1,animate=T)
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

  #summary crossICC result----
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

  #Plot functions----
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
        sih<-fuck[[input$iterslided]]$clusters$silhouette
        max.sliw<-which.max(max(sih[,3])) + 1
        colorlength <- 3
        if(length(unique(sih[,1]))>3){
          colorlength <- length(unique(sih[,1]))
        }
        color.list<-brewer.pal(colorlength, "Set2")
        plot(sih,col=color.list[1:max.sliw])
      })
      output$clusterexpress<-renderPlot({
        validate(
          need(!is.null(InterationResult()), "Please upload a correct CrossICC output file in RDA format, which can be found at default output path of CrossICC function or user defined path.")
        )
        fuck<-InterationResult()
        grid::grid.newpage()
        grid::grid.draw(fuck[[input$iterslided]]$heatmaps[[input$SelectPL]]$gtable)
      })
      output$ssGSEAmatrix<-renderDataTable({
        validate(
          need(!is.null(InterationResult()), "Please upload a correct CrossICC output file in RDA format, which can be found at default output path of CrossICC function or user defined path.")
        )
        fuck<-InterationResult()
        ssGSEA.list<-ssGSEA(fuck[[input$SelectPL]], fuck[[input$iterslided]]$gene.signature, CrossICC.object[[input$iterslided]]$unioned.genesets)
        ssGSEA.list[[1]]
      })
      # output$ssGSEAheatmap-renderPlot({
      #   validate(
      #     need(!is.null(InterationResult()), "Please upload a correct CrossICC output file in RDA format, which can be found at default output path of CrossICC function or user defined path.")
      #   )
      #   fuck<-InterationResult()
      #   ssGSEA.list<-ssGSEA(fuck[[input$SelectPL]], fuck[[input$iterslided]]$gene.signature, CrossICC.object[[input$iterslided]]$unioned.genesets)
      #   ssGSEA.list[[2]]
      # })
  #Download functions----
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
          sih<-fuck[[input$iterslided]]$clusters$silhouette
          max.sliw<-which.max(max(sih[,3])) + 1
          color.list<-brewer.pal(length(unique(sih[,1])), "Set2")
          plot(sih,col=color.list[1:max.sliw])
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
# Predict panel functions ----

      predict.inputdata<- reactive({
        inFile <- input$file2
        data<-NULL
        if (!is.null(inFile)){
          data<-read.csv(data,header=T,row.names=1,check.names = F)
        }
        data
      })
      output$predictInputDataSummary<-renderDataTable({
        predict.inputdata()
      })
})
