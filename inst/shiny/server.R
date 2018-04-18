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
    if(is.null(InterationResult())){
      tagList(div())
    }else{
      df=InterationResult()
      platformnamelist<-names(df[[input$iterslided]]$platforms)
      tagList(
        selectInput("SelectPL", "SelectPlatform", choices=platformnamelist, selected = platformnamelist[1], multiple = FALSE)
      )
    }
  })
  #crossICC outputData ----
      inputdata<-  reactive({
        outputdata<- readRDS(file = path.expand('~/CrossICC.object.rds'))
        inFile <- input$file1
        CrossICC.object<-NULL
        if (!is.null(inFile)){
          CrossICC.object<-readRDS(file = inFile$datapath)
        }
        switch(input$dataset,
               "default" = outputdata,
               "upload" = CrossICC.object
        )
      })
  #Interation CrossICC
      InterationResult <- eventReactive(input$submit, {
        CrossICC.object=inputdata()
        CrossICC.object
      })

        summary.data<-eventReactive(input$submit,{
          temp.summary <-  CrossICC::summary.CrossICC(inputdata(),iteration = input$iterslided)
          return(temp.summary)
        })

  #Summary crossICC result----
      output$OutputResultSignature <- renderDataTable({

        temp.summary <-  CrossICC::summary.CrossICC(InterationResult(),iteration = input$iterslided)
        df<-temp.summary$gene.signatures
        colnames(df)<-c("Cluster","Feature")
        df
        })
      output$OutputClusterResult <- downloadHandler(
        filename = function() {
          paste("Final_cluster_result", Sys.time(), '.csv', sep='')
        },
        content = function(file) {
          temp.summary <-  CrossICC::summary.CrossICC(InterationResult(),iteration = input$iterslided)

          write.csv(temp.summary$cluster, file)

        },
        contentType = 'text/csv'

      )

  # Render arguments matrix----
      output$outputArguments <- renderTable({
        validate(
          need(!is.null(InterationResult()), "Press submit button")
        )
        crossICC.object<-InterationResult()
        validate(
          need(!is.null(input$iterslided), "Press submit button")
        )

        arg.list.2<-crossICC.object[[input$iterslided]]$arg.list
        tempname<-names(arg.list.2)
        tempname[1]="Input"
        df<-data.frame(Parameter=tempname,Value=unlist(as.character(arg.list.2), use.names=FALSE))
        df
      })


  #Plot functions----
      output$superclusterPlot<-renderPlot({
        validate(
          need(!is.null(InterationResult()), "Press submit button")
        )
        validate(
          need(!is.null(input$iterslided), "Press submit button")
        )
        crossICC.object<-InterationResult()

        Sys.sleep(1)
        plot_balanced_heatmap(crossICC.object[[input$iterslided]]$clusters$all.k)

      })
      output$Silhouette<-renderPlot({
        Sys.sleep(1)
        validate(
          need(!is.null(InterationResult()), "Press submit button")
        )
        validate(
          need(!is.null(input$iterslided), "Press submit button")
        )
        crossICC.object<-InterationResult()
        sih<-crossICC.object[[input$iterslided]]$clusters$silhouette
        plot_sihouttle_with_crossICCout(sih)

      })
      output$clusterexpress<-renderPlot({
        Sys.sleep(1)
        validate(
          need(!is.null(InterationResult()), "Press submit button")
        )
        validate(
          need(!is.null(input$iterslided), "Press submit button")
        )
        validate(
          need(!is.null(input$SelectPL), "Press submit button")
        )
        crossICC.object<-InterationResult()
        #plot heatmap
        # get data
        plot.matrix<-as.data.frame(crossICC.object[[input$iterslided]]$platforms[[input$SelectPL]])
        platform.names <- names(crossICC.object[[input$iterslided]]$platforms)
        index <- which(platform.names %in% input$SelectPL)
        cluster.table<-crossICC.object[[input$iterslided]]$clusters$clusters
        gsig<-crossICC.object[[input$iterslided]]$gene.order
        #plot
        plot_expression_heatmap_with_cluster(plot.matrix,cluster.table,gsig)

      })
      output$ssGSEAmatrix<-renderDataTable({
        validate(
          need(!is.null(InterationResult()), "Press submit button")
        )
        validate(
          need(!is.null(input$iterslided), "Press submit button")
        )
        crossICC.object<-InterationResult()
        ssGSEA.list<-ssGSEA(crossICC.object[[input$iterslided]]$platforms[[input$SelectPL]], crossICC.object[[input$iterslided]]$gene.signature, crossICC.object[[input$iterslided]]$unioned.genesets,cluster = crossICC.object[[input$iterslided]]$clusters$clusters)
        ssGSEA.list[[1]]
      })

  #Download functions----
      output$DownloadSuperclusterPlot<-downloadHandler(
        filename = function() {
          paste("SuperclusterPlot_", Sys.time(), '.pdf', sep='')
        },
        content = function(file) {
          pdf(file)
          crossICC.object<-InterationResult()
          plot_balanced_heatmap(crossICC.object[[input$iterslided]]$clusters$all.k)
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
          crossICC.object<-InterationResult()
          sih<-crossICC.object[[input$iterslided]]$clusters$silhouette
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
          crossICC.object<-InterationResult()
          plot.matrix<-as.data.frame(crossICC.object[[input$iterslided]]$platforms[[input$SelectPL]])
          cluster.table<-crossICC.object[[input$iterslided]]$clusters$clusters
          gsig<-crossICC.object[[1]]$gene.order
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
          crossICC.object<-InterationResult()
          plot.matrix<-as.data.frame(crossICC.object[[input$iterslided]]$platforms[[input$SelectPL]])
          write.csv(plot.matrix, file)

        },
        contentType = 'text/csv'
      )
      output$geneSignature<-downloadHandler(
        filename = function() {
          paste("GeneSignarure", Sys.time(), '.csv', sep='')
        },
        content = function(file) {
          temp.summary <-  CrossICC::summary.CrossICC(InterationResult(),iteration = input$iterslided)
          df<-temp.summary$gene.signatures
          colnames(df)<-c("Cluster","Feature")
          df
          write.csv(df, file)

        },
        contentType = 'text/csv'
      )
# Predict panel functions ----

      predict.inputdata<- reactive({
        inFile <- input$file2
        data<-NULL
        if (!is.null(inFile)){
          data<-read.csv(inFile$datapath,header=T,row.names=1,check.names = F)
        }
        data
      })
      output$predictInputData<-renderDataTable({
        predict.inputdata()
      })
      get_predict_result<-reactive({
        validate(
          need(!is.null(InterationResult()), "Press submit button")
        )
        validate(
          need(!is.null(predict.inputdata()), "Press upload data to predict")
        )
        if(input$submit2==0){
        return(NULL)
        }
        predict.data<-predict.inputdata()
        crossICC.object<-InterationResult()
        #validation.Data shoud be format features in rows and samples in columns
        maxIter<-length(crossICC.object)
        crossICC.object.summary<-summary.CrossICC(crossICC.object)
        # get centroid
        train.centroid<-cluster.centroid(crossICC.object[[maxIter]]$platforms[[1]],crossICC.object[[maxIter]]$gene.signature,crossICC.object.summary$clusters)
        #prediction
        vali.predict.bycentroid<-centroid2exp(train.centroid,predict.data)
        #get prediction result
        vali.predict.bycentroid.cluter<-vali.predict.bycentroid$cluster
        vali.predict.normalized.matrix<-vali.predict.bycentroid$normalized.matrix
        return(list(cluster=vali.predict.bycentroid.cluter,matrix=vali.predict.normalized.matrix))
      })
      # predict heatmap for replication
      output$predictHeatmap<-renderPlot({

        predict.list<-get_predict_result()
        cfDNA.crossICC<-InterationResult()
        max.iter<-length(cfDNA.crossICC)
        cluster<-predict.list$cluster
        plot.matrix<-predict.list$matrix
        plot_expression_heatmap_with_cluster(plot.matrix,cluster,cfDNA.crossICC[[max.iter]]$gene.order)

      })
      output$DownloadPredictHeatmap<-downloadHandler(
        filename = function() {
          paste("PredictHeatmap_", Sys.time(), '.pdf', sep='')
        },
        content = function(file) {
          pdf(file)
          predict.list<-get_predict_result()
          cfDNA.crossICC<-InterationResult()
          max.iter<-length(cfDNA.crossICC)
          cluster<-predict.list$cluster
          plot.matrix<-predict.list$matrix
          validate(
            need(!is.null(plot.matrix), "Press the prediction button")
          )
          plot_expression_heatmap_with_cluster(plot.matrix,cluster,cfDNA.crossICC[[max.iter]]$gene.order)
          dev.off()
        },
        contentType = 'image/pdf'
      )
# clinical correlation analysis
      clinicalRelatedData<-  reactive({
        example<- data()
        inFile <- input$file3
        clinical.df<-NULL
        if (!is.null(inFile)){
          clinical.df<-read.csv(inFile,header=T,check.names = F)
        }
        switch(input$dataset,
               "Default" = outputdata,
               "Upload" = clinical.df
        )
      })
})
