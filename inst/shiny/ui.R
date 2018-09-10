options(spinner.type=4)

#required packages
suppressMessages(library(shiny))
suppressMessages(library(shinydashboard))
suppressMessages(library(shinyWidgets))
suppressMessages(library(shinycssloaders))
suppressMessages(library(DT))


#sider bar----

sidebar <- dashboardSidebar(width = 300,
                            sidebarMenu(id="sidername",
                                        menuItem("Home", tabName = "home", icon = icon("home")),
                                        menuItem("CrossICC result viewer", icon = icon("eye"),tabName = "CrossICC"),
                                        menuItem("Allocate New Sample", icon = icon("location-arrow"),tabName = "predict" ),
                                        menuItem('Correlation analysis', tabName = 'correlation',icon = icon("th",lib = "glyphicon")),
                                        menuItem('ssGSEA analysis', tabName = 'ssgsea',icon = icon("bar-chart")),
                                        menuItem('Survival analysis', tabName = 'Survival',icon = icon("line-chart")),
                                        menuItem("Help", tabName = "Help", icon = icon("question"))
                            )
)


#bodyHome ----

bodyHome <- tabItem("home",
                    fluidRow(
                      box(
                        width = 12,
                        status = "success",
                        solidHeader = TRUE,
                        title = strong("Wellcome to the CrossICC reporter"),
                        includeMarkdown("Readmeshiny.md")
                      )
                    ),

                    fluidRow(
                      box(
                        width = 4,
                        status = "info",
                        solidHeader = FALSE,
                        title =  strong("CrossICC result viewer"),
                        tags$img(src = "images/workflow.png")
                      ),
                      box(
                        width = 4,
                        status = "info",
                        solidHeader = FALSE,
                        title =  strong("Predicting new samples")
                      ),
                      box(
                        width = 4,
                        status = "info",
                        solidHeader = FALSE,
                        title =  strong("Cancer related analysis functions")
                      )
                    ),
                    fluidRow(
                      box(
                        width = 12,
                        status = "warning",
                        solidHeader = TRUE,
                        title =  strong("Citation"),
                        h4("If you use CrossICC, please cite this paper:"),
                        p("Preparing in process")
                      )
                    ))

#bodyCrossICC ----

bodyCrossICC <- tabItem("CrossICC",
                        fluidRow(
                          #setting panel
                          box(
                            title = div(shiny::icon("gear"), "Data upload & Configuration", inline =T),
                            width = 4,
                            tabBox (
                              height = "100%", width = "100%",
                              id = "crossICCsettingPanel",
                              title = "",
                              side = "left",
                              selected = "crs01",
                              tabPanel(
                                title = div(icon("book"), "Upload"),
                                value = "crs01",
                                  radioButtons(
                                    "dataset",
                                    strong("Loading CrossICC output"),
                                    inline = TRUE,
                                    c(Default = "default", Upload = "upload"),
                                    selected = 'default'
                                  ),
                                  conditionalPanel(
                                    condition = "input.dataset == 'upload'",
                                    fileInput(
                                      'file1','CrossICC output data in RDS format',accept = c('application/rds', '.rds')
                                    )
                                  ),
                                  actionBttn("submit", div(
                                    "Click ME to visualize result ",
                                    icon("hand-right", lib = "glyphicon")
                                  )),

                                  div(align = 'center',h4("Your running parameters "), tableOutput("outputArguments"))
                                  ),
                              tabPanel(
                                title = div(icon("book"), "Setting"),
                                value = "crs02",
                                fluidRow(
                                  switchInput(
                                    inputId = "clusterRow",
                                    label = "Cluster by Rows",
                                    value = FALSE
                                  ),
                                  switchInput(
                                    inputId = "showRowNames",
                                    label = "show row names",
                                    value = FALSE
                                  )
                                ),
                                sliderInput("cross_size", "Zoom in/out graph:",min = 500, max = 1000, value = 800)
                              )
                            )
                          ),
                          box(
                            width = 8,
                            tabBox (
                              height = "100%", width = "100%",
                              id = "crossICCresultPanel",
                              title = "",
                              side = "right",
                              selected = "cr01",
                              tabPanel(
                                title = div(icon("book"), "Summary"),
                                value = "cr01",
                                dataTableOutput("OutputResultSignature"),
                                box(
                                  width = NULL,  status = "success",
                                  h4("Dowload Clustering result"),
                                  downloadButton("OutputClusterResult"),
                                  h4("Gene signature for each cluster"),
                                  downloadButton("geneSignature", label = "Download Matrix")
                                  )

                              ),
                              tabPanel(
                                title = div(icon("th", lib = "glyphicon"), "Super Clustering"),
                                value = "cr02",

                                withSpinner(plotOutput("superclusterPlot", height = "100%"), color = "black"),
                                box(
                                  width = NULL,  status = "success",
                                  checkboxGroupInput("DownloadSuperclusterPlot_check","Choose file type to download:",
                                    c("png", "pdf", "tiff"),inline = TRUE
                                  ),
                                  downloadButton('DownloadSuperclusterPlot', 'Download')
                                )
                              ),
                              tabPanel(
                                title = div(icon("signal", lib = "glyphicon"), "Silhouette Result"),
                                value = "cr03",

                                withSpinner(plotOutput("Silhouette",height = "100%"), color = "black")
                                ,
                                box(
                                  width = NULL,  status = "success",
                                  checkboxGroupInput("DownloadSilhouette_check","Choose file type to download:",
                                                     c("png", "pdf", "tiff"),inline = TRUE
                                  ),
                                  downloadLink('DownloadSilhouette', 'Download')
                                  )

                              ),
                              tabPanel(
                                title = div(icon("book"), "Expression heatmap by signagure"),
                                value = "cr04",


                                uiOutput("expressionHeatmapSelectPlatform"),
                                withSpinner(plotOutput("clusterexpress",height = "100%"), color = "black"),
                                box(
                                  width = NULL,  status = "success",
                                  checkboxGroupInput("DownloadClusterexpressPlot_check","Choose file type to download:",
                                                     c("png", "pdf", "tiff"),inline = TRUE
                                  ),
                                  downloadLink('DownloadClusterexpressPlot', 'Download'),
                                  downloadButton('DownloadClusterExpressMatrix', 'Download Matrix')
                                )
                              ),
                              tabPanel(
                                title = div(icon("book"), "Iteration Record"),
                                value = "cr05",
                                plotOutput("IterationPlot", height = "800px"),
                                box(
                                  width = NULL,  status = "success",
                                  checkboxGroupInput("DownloadIterationPlot_check","Choose file type to download:",
                                                     c("png", "pdf", "tiff"),inline = TRUE
                                  )
                                # to be added
                                )
                              )
                            )
                          )
                        ))

#bodyPredict----

bodyPredict <- tabItem(
  "predict",
  h2("Allocate new samples with their expression value"),
  fluidRow(
    #setting panel
    box(
      title = "Input your data set for prediction",
      background = "black",
      width = 4,
      radioButtons(
        "dataset2",
        strong("Data to predict"),
        inline = TRUE,
        c(Default = "default", Upload = "upload"),
        selected = 'default'
      ),
      conditionalPanel(
        condition = "input.dataset2 == 'upload'",
        fileInput(
          'file2',
          'Input dataset in matrix file',
          accept = c('text/txt', '.rds')
        )
      ),
      actionBttn("submit2", "Submit")
    ),
    tabBox (
      id = "PredictResultPanel",
      title = h3("Analysis"),
      width = 8,
      side = "right",
      selected = "pre01",
      tabPanel(
        title = div(icon("book"), "PredictResult"),
        value = "pre01",

        downloadLink('DownloadPredictHeatmap', 'Download PDF'),
        withSpinner(plotOutput("predictHeatmap", height = "800px"), color = "black")

      )
    )
  )

)

# bodyCorrelation----

bodyCorrelation <- tabItem(
  "correlation",
  h2("Correlation analysis of the two cluster system"),
  fluidRow(
    valueBoxOutput("getRAbox", width = 4),
    valueBoxOutput("getARIbox", width = 4),
    valueBoxOutput("getJaccarddox", width = 4)

  ),
  fluidRow(
    tabBox (
      width = 4,
      selected = "caInput01",
      side = "left",
      tabPanel(
        title = div(icon("table"), "Data Input"),
        value = "caInput01",
        radioGroupButtons(
          "data3",
          label = strong("clinical Dataset"),
          choices = c(Default = "Default", Upload = "Upload"),
          selected = 'Default',
          status = "primary"
        ),
        conditionalPanel(
          condition = "input.data3 == 'Upload'",
          fileInput(
            'file3',
            'Input dataset in matrix file',
            accept =
              c('text/txt', '.rds')
          )
        ),
        uiOutput("VariableSelectionUI1"),
        uiOutput("VariableSelectionUI2"),
        actionBttn("submit3", "Submit", style = "unite")

      )
    ),
    tabBox (
      width = 8,
      id = "clinicalResultPanel",
      title = h3("Analysis"),
      side = "right",
      selected = "ca02",
      tabPanel(
        title = div(icon("book"), "Read Me"),
        value = "ca01",
        p("Write introduction here")

      ),
      tabPanel(
        title = div(icon("book"), "Data"),
        value = "ca02",
        h3("Input Data "),
        dataTableOutput("summaryCorrelationData")
      ),
      tabPanel(title = div(icon("book"), "Result"), value =
                 "ca03")
      ,
      tabPanel(
        title = div(icon("book"), "Overlap Analysis"),
        value = "ca04",
        h3("Jaccard Index matrix"),
        h4("Overlap")

      )
    )
  )
)

# bodySsGSEA ----
bodySsGSEA <- tabItem("ssgsea",
                      h2("ssGSEA analysis"),
                      fluidRow(
                        #setting panel
                        box(
                          title = "Input dataset",
                          background = "black",
                          width = 4
                        ),
                        tabBox (width = 8)
                      ))

# bodySurvival ----
bodySurival <- tabItem(
  "Survival",
  h2("Survival evaluation of the clustering system"),
  fluidRow(
    #setting panel
    box(
      title = "Input dataset",
      background = "black",
      width = 4
    ),
    tabBox (width = 8)
  )
)

#-----------






#Main function----
shinyUI(dashboardPage(skin = "black",


  dashboardHeader(title = "CrossICC: iterative consensus clustering of cross-platform data",
                  titleWidth = 600
                  ),
  sidebar,
  dashboardBody(
    tags$head(
      tags$style(HTML("
      .shiny-output-error-validation {
                      color: brown;
                      }
                      ")),
      tags$link(rel = "stylesheet", type = "text/css", href = "css/main.css")
    ),
    tabItems(
      bodyHome,
      bodyCrossICC,
      bodyPredict,
      bodyCorrelation,
      bodySsGSEA,
      bodySurival

    )
  )
))
