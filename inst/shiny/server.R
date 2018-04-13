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
source("globalfunctions/plotFunctions.R")
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
      platformnamelist<-names(df[[input$iterslided]]$platforms)
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

  # Render arguments matrix----
      output$outputArguments <- renderTable({
        fuck<-InterationResult()
        dt<-fuck[[input$iterslided]]$arg.table
        dt
      })
  #Plot functions----
      output$superclusterPlot<-renderPlot({
        validate(
          need(!is.null(InterationResult()), "Please upload a correct CrossICC output file in RDA format, which can be found at default output path of CrossICC function or user defined path.")
        )
        fuck<-InterationResult()
        plot_balanced_heatmap(fuck[[input$iterslided]]$clusters$all.k)

      })
      output$Silhouette<-renderPlot({
        validate(
          need(!is.null(InterationResult()), "Please upload a correct CrossICC output file in RDA format, which can be found at default output path of CrossICC function or user defined path.")
        )
        fuck<-InterationResult()
        sih<-fuck[[input$iterslided]]$clusters$silhouette
        plot_sihouttle_with_crossICCout(sih)

      })
      output$clusterexpress<-renderPlot({
        validate(
          need(!is.null(InterationResult()), "Please upload a correct CrossICC output file in RDA format, which can be found at default output path of CrossICC function or user defined path.")
        )
        fuck<-InterationResult()
        #plot heatmap
        # get data
        plot.matrix<-as.data.frame(fuck[[input$iterslided]]$platforms[[input$SelectPL]])
        platform.names <- names(fuck[[input$iterslided]]$platforms)
        index <- which(platform.names %in% input$SelectPL)
        cluster.table<-fuck[[input$iterslided]]$clusters$clusters
        gsig<-fuck[[input$iterslided]]$sorted.gene.list[[index]]
        #plot
        plot_expression_heatmap_with_cluster(plot.matrix,cluster.table,gsig)

      })
      output$ssGSEAmatrix<-renderDataTable({
        validate(
          need(!is.null(InterationResult()), "Please upload a correct CrossICC output file in RDA format, which can be found at default output path of CrossICC function or user defined path.")
        )
        fuck<-InterationResult()
        ssGSEA.list<-ssGSEA(fuck[[input$SelectPL]], fuck[[input$iterslided]]$gene.signature, fuck[[input$iterslided]]$unioned.genesets)
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
          plot_balanced_heatmap(fuck[[input$iterslided]]$clusters$all.k)
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
          plot_sihouttle_with_crossICCout(sih)
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
          plot.matrix<-as.data.frame(fuck[[input$iterslided]]$platforms[[input$SelectPL]])
          cluster.table<-fuck[[input$iterslided]]$clusters$clusters
          gsig<-fuck[[1]]$gene.signature
          #plot
          plot_expression_heatmap_with_cluster(plot.matrix,cluster.table,gsig)
          dev.off()
        },
        contentType = 'image/pdf'
      )
      output$DownloadClusterExpressMatrix<-downloadHandler(
        filename = function() {
          paste("Clusterexpress_", Sys.time(), '.csv', sep='')
        },
        content = function(file) {
          pdf(file)
          plot.matrix<-as.data.frame(fuck[[input$iterslided]]$platforms[[input$SelectPL]])
          write.csv(plot.matrix, file)

        },
        contentType = 'text/csv'
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
