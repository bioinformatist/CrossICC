options(spinner.type=4)

#required packages
suppressMessages(library(shiny))
suppressMessages(library(shinydashboard))
suppressMessages(library(shinyWidgets))
suppressMessages(library(shinycssloaders))
suppressMessages(library(DT))

# library(surival)d
# suppressMessages(library(shinyIncubator))
# suppressMessages(library(plotly))
# suppressMessages(library(heatmaply))

#Main function
shinyUI(dashboardPage(skin = "black",


  dashboardHeader(title = "CrossICC: iterative consensus clustering of cross-platform gene expression data",
                  titleWidth = 600
                  ),
  # sider bar ----
  dashboardSidebar(sidebarMenu(
    menuItem("Home", tabName = "home", icon = icon("home")),
    menuItem("Analysis", icon = icon("search"),tabName = "analysis",
             collapsible = T,
             menuSubItem('CrossICC result', tabName = 'CrossICC'),
             menuSubItem('Predict New Sample', tabName = 'predict'),
             menuSubItem('Correlation analysis', tabName = 'correlation')
    ),
    menuItem("Help", tabName = "Help", icon = icon("question"))
  )),
  #body elements ----
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
      # home page----
      tabItem("home",
              fluidRow(
                box(
                  width = 6, status = "success", solidHeader = TRUE,collapsible = TRUE,
                  title = "Read Me",
                  includeMarkdown("Readmeshiny.md")
                ),
                box(
                  width = 6, status = "info", solidHeader = TRUE,collapsible = TRUE,
                  title = "Workflow",
                  tags$img(src="images/workflow.png")
                )
              )
      ),

      # inout panel ----
      tabItem("CrossICC",
              fluidRow(
                #setting panel
                  box(
                   title = div(shiny::icon("gear"),"Control Panel for Example data", inline=T),width = 4, background = "black",
                   radioButtons(
                     "dataset",
                     strong("Loading CrossICC output"),inline=TRUE,
                     c(Default = "default", Upload = "upload"),
                     selected = 'default'
                   ),
                   conditionalPanel(condition = "input.dataset == 'upload'",
                                    fileInput('file1', 'CrossICC output data in RDS format',
                                              accept=c('application/rds', '.rds'))
                   ),
                   actionBttn("submit",div("Click ME to visualize result ",icon("hand-right",lib = "glyphicon")),style = "float"),
                   h4("Your running parameters "),
                   div(align='center',tableOutput("outputArguments"))
                  ),
                  tabBox (id="crossICCresultPanel",title="", width = 8,side = "right",
                          selected = "cr01",
                          tabPanel(title=div(icon("book"),"Summary"),value="cr01",
                                   h4("Dowload Clustering result"),
                                   downloadButton("OutputClusterResult"),
                                   h4("Gene signature for each cluster"),
                                   downloadButton("geneSignature",label = "Download Matrix"),
                                   dataTableOutput("OutputResultSignature")

                          ),
                          tabPanel(title=div(icon("th",lib = "glyphicon"),"Super Clustering"),value="cr02",
                                   downloadLink('DownloadSuperclusterPlot', 'Download PDF'),
                                   withSpinner(plotOutput("superclusterPlot",height = "800px"),color = "black")
                          ),
                          tabPanel(title=div(icon("signal",lib = "glyphicon"),"Silhouette Result"),value="cr03",
                                   downloadLink('DownloadSilhouette', 'Download PDF'),
                                   withSpinner( plotOutput("Silhouette",height = "800px"),color = "black")
                          ),
                          tabPanel(title=div(icon("book"),"Expression heatmap by signagure"),value="cr04",
                                   downloadLink('DownloadClusterexpressPlot', 'Download PDF'),
                                   downloadButton('DownloadClusterExpressMatrix', 'Download Matrix'),
                                   materialSwitch(inputId = "clusterRow",
                                                  label = "Cluster by Rows", status = "primary",
                                                  right = FALSE),
                                   materialSwitch(inputId = "showRowNames",
                                                  label = "show row names", status = "primary",
                                                  right = FALSE),
                                   uiOutput("expressionHeatmapSelectPlatform"),
                                   withSpinner(plotOutput("clusterexpress",height = "800px"),color = "black")
                          ),
                          tabPanel(title=div(icon("book"),"Iteration Record"),value="cr05",
                                   plotOutput("IterationPlot",height = "800px")
                          )
                  )
            )
      ),
      tabItem("predict",
              fluidRow(
                #setting panel
                box(
                  title = "Input your data set for prediction", background = "black",width = 4,
                  radioButtons(
                    "dataset2",
                    strong("Data to predict"),inline=TRUE,
                    c(Default = "default", Upload = "upload"),
                    selected = 'default'
                  ),
                  conditionalPanel(condition = "input.dataset2 == 'upload'",
                                   fileInput('file2', 'Input dataset in matrix file',
                                             accept=c('text/txt', '.rds'))
                  ),
                  actionBttn("submit2","Submit")
                ),
                tabBox (id="PredictResultPanel",title=h3("Analysis"),width = 8, side = "right",
                        selected = "pre01",
                        tabPanel(title=div(icon("book"),"PredictResult"),value="pre01",

                                 downloadLink('DownloadPredictHeatmap', 'Download PDF'),
                                 withSpinner(plotOutput("predictHeatmap",height = "800px"),color = "black")

                        )
                )
              )

              ),
      tabItem("correlation",
              fluidRow(
                #setting panel
                box(
                  title = "Phenotype Data input",status = "success",width = 4, background = "black",
                  radioGroupButtons(
                    "data3",
                    label = strong("clinical Dataset"),
                    choices = c(Default = "Default", Upload = "Upload"),
                    selected = 'Default',status = "primary"
                  ),
                  conditionalPanel(condition = "input.data3 == 'Upload'",
                                   fileInput('file3', 'Input dataset in matrix file',
                                             accept=c('text/txt', '.rds'))
                  ),
                  uiOutput("VariableSelectionUI1"),
                  uiOutput("VariableSelectionUI2"),
                  actionBttn("submit3","Submit", style = "unite")
                ),


                       valueBoxOutput("getRAbox",width = 2),
                       valueBoxOutput("getARIbox",width = 2)

                       # tabBox (id="clinicalResultPanel",title=h3("Analysis"), side = "right",
                       #         selected = "ca02",
                       #         tabPanel(title=div(icon("book"),"Read Me"),value="ca01",
                       #                  p("Write introduction here")
                       #
                       #         ),
                       #         tabPanel(title=div(icon("book"),"Data"),value="ca02",
                       #                  h3("Input Data "),
                       #                  dataTableOutput("summaryCorrelationData")
                       #         ),
                       #         tabPanel(title=div(icon("book"),"Result"),value="ca03"
                       #
                       #
                       #
                       #         )
                       #         ,
                       #         tabPanel(title=div(icon("book"),"Overlap Analysis"),value="ca04",
                       #                  h3("Jaccard Index matrix"),
                       #                  h4("Overlap")
                       #
                       #         )
                       # )


              )

              )

      # analysis  panel ----
    )
    # Boxes need to be put in a row (or column)



  )
))
