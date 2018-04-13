

#required packages
suppressMessages(library(shiny))
suppressMessages(library(shinydashboard))
suppressMessages(library(shinyWidgets))
suppressMessages(library(DT))
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
                   title = div(shiny::icon("gear"),"Control Panel for Example data "),width = 4, background = "black",
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
                   actionButton("submit","Click ME to visualize result "),
                   uiOutput("interationNumberForplot"),
                   tableOutput("outputArguments")
                  ),
                  tabBox (id="crossICCresultPanel",title=div(icon("hand-right",lib = "glyphicon"),h3("Data Exploration")), width = 8,side = "right",
                          selected = "cr01",
                          tabPanel(title=div(icon("book"),"Summary"),value="cr01",
                                   h3("Sample clusters"),
                                   dropdownButton(
                                     circle = TRUE, status = "primary", icon = icon("gear"), width = "300px",
                                     tooltip = tooltipOptions(title = "Click to see inputs !")
                                   ),
                                   verbatimTextOutput("OutputClusterResult"),
                                   h3("Gene signature for each cluster"),
                                   dataTableOutput("OutputResultSignature")
                          ),
                          tabPanel(title=div(icon("th",lib = "glyphicon"),"Super Clustering"),value="cr02",
                                   downloadLink('DownloadSuperclusterPlot', 'Download PDF'),
                                   plotOutput("superclusterPlot",height = "800px")
                          ),
                          tabPanel(title=div(icon("signal",lib = "glyphicon"),"Silhouette Result"),value="cr03",
                                   downloadLink('DownloadSilhouette', 'Download PDF'),
                                   plotOutput("Silhouette",height = "800px")
                          ),
                          tabPanel(title=div(icon("book"),"Expression heatmap by signagure"),value="cr04",
                                   downloadLink('DownloadClusterexpressPlot', 'Download PDF'),
                                   downloadLink('DownloadClusterExpressMatrix', 'Download Matrix'),
                                   uiOutput("expressionHeatmapSelectPlatform"),
                                   plotOutput("clusterexpress",height = "800px")
                          ),
                          tabPanel(title=div(icon("book"),"ssGSEA"),value="cr05",
                                   plotOutput("ssGSEAheatmap",height = "800px"),
                                   dataTableOutput("ssGSEAmatrix")
                          )
                  )
            )
      ),
      tabItem("predict",
              fluidRow(
                #setting panel
                box(
                  title = "Input your data set for prediction",solidHeader = TRUE,status = "success",width = 4,
                  radioButtons(
                    "dataset2",
                    strong("Mutation Dataset"),inline=TRUE,
                    c(Default = "default", Upload = "upload"),
                    selected = 'default'
                  ),
                  conditionalPanel(condition = "input.dataset2 == 'upload'",
                                   fileInput('file2', 'Input dataset in matrix file',
                                             accept=c('text/txt', '.rds'))
                  ),
                  actionButton("submit2","Submit")
                ),
                tabBox (id="PredictResultPanel",title=h3("Analysis"),width = 8, side = "right",
                        selected = "pre01",
                        tabPanel(title=div(icon("book"),"Read Me"),value="pre01",
                                 p("Write introduction here")
                        ),
                        tabPanel(title=div(icon("book"),"Data Input"),value="pre02",
                                 h3("Summary of input dataset "),
                                 dataTableOutput("predictInputDataSummary")

                        ),
                        tabPanel(title=div(icon("book"),"Predict Result"),value="pre03"

                        )
                )
              )

              ),
      tabItem("correlation")

      # analysis  panel ----
    )
    # Boxes need to be put in a row (or column)



  )
))
