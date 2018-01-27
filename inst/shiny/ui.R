

#required packages
suppressMessages(library(shiny))
suppressMessages(library(shinydashboard))

suppressMessages(library(DT))
# suppressMessages(library(shinyIncubator))
# suppressMessages(library(plotly))
# suppressMessages(library(heatmaply))

#Main function
shinyUI(dashboardPage(skin = "black",
  dashboardHeader(title = "CrossICC: iterative consensus clustering of cross-platform gene expression data",
                  titleWidth = 600),
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
    tabItems(
      # home page----
      tabItem("home",
              box(
                width = 8, status = "info", solidHeader = TRUE,collapsible = TRUE,
                title = "Workflow",
                tags$img(src="images/workflow.png")
              )

              ),

      # inout panel ----
      tabItem("CrossICC",
              fluidRow(
                  box(
                   title = "Control Panel for Example data ",solidHeader = TRUE,status = "success",width = 4,
                   numericInput("MaxInterNum","Max iterater number",value=1000,min=100,max=1000,step=100),
                   actionButton("submit","Submit"),
                   h3("Run result showing here"),
                    uiOutput("interationNumberForplot")
                  ),
                  tabBox (id="crossICCresultPanel",title=h3("Data Exploration"),width = 8, side = "right",
                          selected = "cr01",
                          tabPanel(title=div(icon("book"),"Super Clustering"),value="cr01",
                                   plotOutput("superclusterPlot")
                          ),
                          tabPanel(title=div(icon("book"),"Silhouette Result"),value="cr02",
                                   plotOutput("Silhouette")
                          ),
                          tabPanel(title=div(icon("book"),"Expression heatmap by signagure"),value="cr03",
                                   uiOutput("expressionHeatmapSelectPlatform"),
                                   plotOutput("clusterexpress")
                          )
                  )
                )
              ),
      tabItem("predict"),
      tabItem("correlation")

      # analysis  panel ----
    )
    # Boxes need to be put in a row (or column)



  )
))
