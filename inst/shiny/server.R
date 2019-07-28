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
suppressMessages(library(reshape2))
suppressMessages(library(ggsci))
suppressMessages(library(ggthemes))
suppressMessages(library(tibble))
source("globalfunctions/plotFunctions.R")
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

#=============================#
# CrossICC panel functions ----
#=============================#

  cross_size <- reactive({
    return(input$cross_size)
  })

  corre_size <- reactive({
    return(input$corre_size)
  })

  predict_size <- reactive({
    return(input$predict_size)
  })

  ssgsea_size <- reactive({
    return(input$ssgsea_size)
  })
  survival_size <- reactive({
    return(input$survival_size)
  })

  #animate bar ---

  #heatmap control option of platform selection ui----
  output$expressionHeatmapSelectPlatform <- renderUI({
    if(is.null(InterationResult())){
      tagList(div())
    }else{
      CrossICC.object=InterationResult()
      platformnamelist<-names(CrossICC.object$platforms)
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
  #Interation CrossICC----
      InterationResult <- eventReactive(input$submit, {
        CrossICC.object=inputdata()
        CrossICC.object
      })

        summary.data<-eventReactive(input$submit,{
          temp.summary <-  CrossICC::summaryCrossICC(inputdata())
          return(temp.summary)
        })

  #Summary crossICC result----
      output$OutputResultSignature <- renderDataTable({

        temp.summary <-  CrossICC::summaryCrossICC(InterationResult())
        df<-temp.summary$gene.signatures
        colnames(df)<-c("Cluster","Feature")
        df
        })
      output$OutputClusterResult <- downloadHandler(
        filename = function() {
          paste("Final_cluster_result", Sys.time(), '.csv', sep='')
        },
        content = function(file) {
          temp.summary <-  CrossICC::summaryCrossICC(InterationResult())

          write.csv(temp.summary$clusters, file)

        },
        contentType = 'text/csv'

      )

  # Render arguments matrix----
      output$outputArguments <- renderTable({
        validate(
          need(!is.null(InterationResult()), "Press submit button")
        )
        crossICC.object<-InterationResult()

        arg.list.2<-crossICC.object$arg.list
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

        crossICC.object<-InterationResult()

        xx<-plot_balanced_heatmap(crossICC.object$clusters$all.k)
        grid::grid.draw(xx$gtable)

      },
        width = cross_size,
        height = cross_size,
        outputArgs = list()
      )
      output$Silhouette<-renderPlot({
        validate(
          need(!is.null(InterationResult()), "Press submit button")
        )

        crossICC.object<-InterationResult()
        sih<-crossICC.object$clusters$silhouette
        plot_sihouttle_with_crossICCout(sih)

      },
        width = cross_size,
        height = cross_size,
        outputArgs = list()
      )

      platform<-reactive({
        input$SelectPL
      })

      getClusterexpress<-reactive({
        validate(
          need(!is.null(InterationResult()), "Press submit button")
        )

        validate(
          need(!is.null(input$SelectPL), "Press submit button")
        )
        crossICC.object<-InterationResult()
        platform.names <- names(crossICC.object$platforms)
        index <- which(platform.names %in% platform())
        plot.matrix<-as.data.frame(crossICC.object$platforms[[index]])
        if(is(crossICC.object$clusters$clusters, "list")){
          cluster.table<-crossICC.object$clusters$clusters[[index]]
        }else{
          cluster.table<-crossICC.object$clusters$clusters
        }
        gsig<-crossICC.object$gene.order[[index]]
        #plot
        xx<-plot_expression_heatmap_with_cluster(plot.matrix,cluster.table,gsig,cluster_row = input$clusterRow,showRowname = input$showRowNames)
        return(xx)
      })
      output$clusterexpress<-renderPlot({

        getClusterexpress()

      },
        width = cross_size,
        height = cross_size,
        outputArgs = list()
      )


      output$IterationPlot<-renderPlot({
        validate(
          need(!is.null(InterationResult()), "Press submit button")
        )
        crossICC.object<-InterationResult()
        plot_iter_with_crossICC(crossICC.object)
      },
        width = cross_size,
        height = cross_size,
        outputArgs = list()
      )


  #Download functions----
      #plot
      output$DownloadSuperclusterPlot<-downloadHandler(

        filename = function() {
          paste("SuperclusterPlot_", Sys.time(), '.',input$DownloadSuperclusterPlotCheck, sep='')
        },
        content = function(file) {
          pixelratio <- 2

          if (input$DownloadSuperclusterPlotCheck == "png")
            png(file, res=72*pixelratio, units = "px")
          else if (input$DownloadSuperclusterPlotCheck == "pdf")
            pdf(file)
          else
            tiff(file)

          crossICC.object<-InterationResult()
          plot_balanced_heatmap(crossICC.object$clusters$all.k)
          dev.off()
        },
        contentType = paste('image/',input$DownloadSuperclusterPlotCheck,sep="")
      )
      output$DownloadSilhouette<-downloadHandler(
        filename = function() {
          paste("Silhouette_", Sys.time(),  '.',input$DownloadSilhouetteCheck, sep='')
        },
        content = function(file) {
          pixelratio <- 2

          if (input$DownloadSilhouetteCheck == "png")
            png(file, res=72*pixelratio, units = "px")
          else if (input$DownloadSilhouetteCheck == "pdf")
            pdf(file)
          else
            tiff(file)
          crossICC.object<-InterationResult()
          sih<-crossICC.object$clusters$silhouette
          plot_sihouttle_with_crossICCout(sih)
          dev.off()
        },
        contentType = paste('image/',input$DownloadSilhouetteCheck,sep="")
      )
      output$DownloadClusterexpressPlot<-downloadHandler(
        filename = function() {
          paste("Clusterexpress_", Sys.time(), '.',input$DownloadClusterexpressPlotCheck, sep='')
        },
        content = function(file) {
          pixelratio <- 2
          if (input$DownloadClusterexpressPlotCheck == "png")
            png(file, res=72*pixelratio, units = "px")
          else if (input$DownloadClusterexpressPlotCheck == "pdf")
            pdf(file)
          else
            tiff(file)

          getClusterexpress()
          dev.off()
        },

        contentType = paste('image/',input$DownloadClusterexpressPlotCheck,sep="")
      )
      #matrix
      output$DownloadClusterExpressMatrix<-downloadHandler(
        filename = function() {
          paste("Clusterexpress_", Sys.time(), '.csv', sep='')
        },
        content = function(file) {
          crossICC.object<-InterationResult()
          plot.matrix<-as.data.frame(crossICC.object$platforms[[input$SelectPL]])
          write.csv(plot.matrix, file)

        },
        contentType = 'text/csv'
      )
      #matrix
      output$geneSignature<-downloadHandler(
        filename = function() {
          paste("GeneSignarure", Sys.time(), '.csv', sep='')
        },
        content = function(file) {
          temp.summary <-  CrossICC::summaryCrossICC(InterationResult())
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
          data<-read.csv(inFile$datapath,header=TRUE,row.names=1,check.names = FALSE)
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
        res.pred<-predictor(predict.data,crossICC.object)
        return(res.pred)
      })
      # predict heatmap for replication
      output$predictHeatmap<-renderPlot({

        predict.list<-get_predict_result()
        cfDNA.crossICC<-InterationResult()
        cluster<-predict.list$cluster
        plot.matrix<-predict.list$matrix
        plot_expression_heatmap_with_cluster(plot.matrix,cluster,cfDNA.crossICC$gene.order)

      })
      output$DownloadPredictHeatmap<-downloadHandler(
        filename = function() {
          paste("PredictHeatmap_", Sys.time(), '.pdf', sep='')
        },
        content = function(file) {
          pdf(file)
          predict.list<-get_predict_result()
          cfDNA.crossICC<-InterationResult()
          cluster<-predict.list$cluster
          plot.matrix<-predict.list$matrix
          validate(
            need(!is.null(plot.matrix), "Press the prediction button")
          )
          plot_expression_heatmap_with_cluster(plot.matrix,cluster,cfDNA.crossICC$gene.order)
          dev.off()
        },
        contentType = 'image/pdf'
      )
# ----------------

# Correlation analysis -----------

      clinicalRelatedData<-  reactive({
        example<-  read.csv(file = path.expand('data/TCGA.COAD.csv'),header = TRUE, row.names = 1)
        inFile <- input$file3
        clinical.df<-NULL
        if (!is.null(inFile)){
          clinical.df<-read.csv(inFile,header=TRUE,check.names = FALSE)
        }
        switch(input$data3,
               "Default" = example,
               "Upload" = clinical.df
        )
      })
      # view data
      output$summaryCorrelationData<-renderDataTable(
        clinicalRelatedData()
      )
      # select data UI
      output$VariableSelectionUI1<-renderUI({
        condition=colnames(clinicalRelatedData())
        conditionVector=as.character(condition)
        selectInput("corAnalysisSelect1", "Variable 1:",choices=conditionVector,selected=conditionVector[1])
      })
      output$VariableSelectionUI2<-renderUI({
        condition=colnames(clinicalRelatedData())
        conditionVector=as.character(condition)
        selectInput("corAnalysisSelect2", "Variable 2:",choices=conditionVector,selected=conditionVector[2])
      })
      # get correlation analysis result
      output$getRAbox<-renderValueBox({
        df<-clinicalRelatedData()
        df<-df[complete.cases(df),]
        x<-input$corAnalysisSelect1
        y<-input$corAnalysisSelect2

        RI<-round(rand.index(df,x,y),digits = 4)
        valueBox(
          "Rand Index",
          RI,
          icon = icon("credit-card")
        )
      })
      output$getARIbox<-renderInfoBox({
        df<-clinicalRelatedData()
        df<-df[complete.cases(df),]
        x<-input$corAnalysisSelect1
        y<-input$corAnalysisSelect2
        ARI<-round(Cal.ARI(df,x,y),digits = 4)

        valueBox(
          "Adjust Rand Index",
          ARI,
          icon = icon("cog",lib = "glyphicon"),
          color = "red"
        )
      })
      output$getJaccarddox<-renderInfoBox({
        df<-clinicalRelatedData()
        df<-df[complete.cases(df),]
        x<-input$corAnalysisSelect1
        y<-input$corAnalysisSelect2
        JI<-round(Cal.ARI(df,x,y),digits = 4)

        valueBox(
          "Jaccard Index",
          JI,
          icon = icon("road",lib = "glyphicon"),
          color = "green"
        )
      })

      getContingencyTable<-reactive({
        df<-clinicalRelatedData()
        df<-df[complete.cases(df),]
        x<-input$corAnalysisSelect1
        y<-input$corAnalysisSelect2
        temp<-data.frame(table(df[,c(x,y)]))
        return(temp)
      })

      output$ContingencyTableRender<-renderDataTable({
        getContingencyTable()
      })

      getcorplot<-reactive({
        df<-clinicalRelatedData()
        df<-df[complete.cases(df),]
        x<-input$corAnalysisSelect1
        y<-input$corAnalysisSelect2
        g<-plotStackBarplot(df,int.vect1 = x,int.vect2 = y,input.theme = input$cor_theme)
        print(g)
      })
      getSankyplot<-reactive({
        df<-clinicalRelatedData()
        df<-df[complete.cases(df),]
        x<-input$corAnalysisSelect1
        y<-input$corAnalysisSelect2
        g<-Sankeyplot(df,int.vect1 = x,int.vect2 = y,input.theme = input$cor_theme)
        print(g)
      })

      #correlation plot
     output$getCorplotRender<-renderPlot({
       getcorplot()
     },
       width = corre_size,
       height = corre_size,
       outputArgs = list()
     )
     output$getSankyPlotRender<-renderPlot({
       getSankyplot()
     },
     width = corre_size,
     height = corre_size,
     outputArgs = list()
     )
      #download plot
     output$DownloadCorrelationPlot<-downloadHandler(
       filename = function() {
         paste("correlattion_", Sys.time(), '.',input$DownloadCorrelationPlot_check, sep='')
       },
       content = function(file) {
         switch (input$DownloadCorrelationPlot_check,
                pdf=pdf(file),
                png=png(file),
                tiff=tiff(file))
         getcorplot()
         validate(
           need(!is.null(plot.matrix), "No data input")
         )
         dev.off()
       },
       contentType =  switch (input$DownloadCorrelationPlot_check,
                              pdf='image/pdf',
                              png='image/pdf',
                              tiff='image/tiff')
     )
     #download plot
     output$DownloadSankeyPlot<-downloadHandler(
       filename = function() {
         paste("sankey_", Sys.time(), '.',input$DownloadCorrelationPlot_check, sep='')
       },
       content = function(file) {
         switch (input$DownloadsankeyPlotCheck,
                 pdf=pdf(file),
                 png=png(file),
                 tiff=tiff(file))
         getSankyplot()
         validate(
           need(!is.null(plot.matrix), "No data input")
         )
         dev.off()
       },
       contentType =  switch (input$DownloadsankeyPlotCheck,
                              pdf='image/pdf',
                              png='image/pdf',
                              tiff='image/tiff')
     )

# -------------------------------

     #GSEA analysis----
     output$ssGSEAmatrix<-renderDataTable({
       validate(
         need(!is.null(InterationResult()), "Press submit button")
       )

       crossICC.object<-InterationResult()
       ssGSEA.list<-ssGSEA(crossICC.object$platforms[[input$SelectPL]], crossICC.object$gene.signature, crossICC.object$unioned.genesets,cluster = crossICC.object$clusters$clusters)
       ssGSEA.list[[1]]
     })
     ssGSEAData<-  reactive({
       example<-  read.csv(file = path.expand('data/surv.test.csv'),header = TRUE,row.names = 1)
       inFile <- input$ssGSEAdatafile
       clinical.df<-NULL
       if (!is.null(inFile)){
         clinical.df<-read.csv(inFile,header=TRUE,check.names = FALSE)
       }
       switch(input$data3,
              "Default" = example,
              "Upload" = clinical.df
       )
     })
##### survival analysis

   SurvivalData<-  reactive({
       example<-  read.csv(file = path.expand('data/surv.test.csv'),header = TRUE,row.names = 1)
       inFile <- input$survivalFile
       clinical.df<-NULL
       if (!is.null(inFile)){
         clinical.df<-read.csv(inFile,header=TRUE,check.names = FALSE)
       }
       switch(input$data3,
              "Default" = example,
              "Upload" = clinical.df
       )
     })
#--------------------------------




})
