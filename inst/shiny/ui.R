

#required packages
suppressMessages(library(shiny))
suppressMessages(library(shiny.semantic))
suppressMessages(library(shinyjs))
suppressMessages(library(shinydashboard))
suppressMessages(library(shinyBS))
suppressMessages(library(DT))
suppressMessages(library(shinyIncubator))
suppressMessages(library(plotly))
suppressMessages(library(heatmaply))

#Main function
shinyUI(dashboardPage(
  dashboardHeader(title = "CrossICC: iterative consensus clustering of cross-platform gene expression data",
                  titleWidth = 600),
  # sider bar ----
  dashboardSidebar(sidebarMenu(
    menuItem("Home", tabName = "home", icon = icon("home")),
    menuItem("Analysis", icon = icon("search"),tabName = "analysis",
             collapsible = T,
             menuSubItem('Input Dataset', tabName = 'input'),
             menuSubItem('Iterater', tabName = 'iter')
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
      tabItem("input",
              fluidRow(
                column(
                  width = 2,
                  #step1----
                  box(
                    title =  div(icon("file-text"),"Example dataset"),solidHeader = TRUE,width = 100,status = "success",
                    radioButtons(
                      "dataset",
                      strong("Mutation Dataset"),inline=T,
                      c(Example = "example"),
                      selected = 'example'
                    ),
                    numericInput("MaxInterNum","Max iterater number",value=1000,min=100,max=1000,step=100),
                    #input manually ----
                    # radioButtons(
                    #   "dataset",
                    #   strong("Mutation Dataset"),inline=T,
                    #   c(Example = "example", Upload = "upload"),
                    #   selected = 'example'
                    # ),
                    # conditionalPanel(condition = "input.dataset == 'upload'",
                    #                  fileInput('file1', 'CSV and Text Document format are supported',
                    #                            accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv'),multiple = T)
                    #                  # ,fileInputSeries()
                    # ),
                    #----
                    actionButton("submit","Submit")
                  ),
                  box(
                    title = "Interation information",solidHeader = TRUE,width = 100,status = "success"

                  ),
                  box(
                    title = "Control panel",solidHeader = TRUE,width = 100,status = "success"
                  )
                ),
                column(
                  width = 4,
                  box(
                    title = "Super Clustering",solidHeader = TRUE,status = "primary"
                  ),
                  box(
                    title = "Silhouette Result",solidHeader = TRUE,status = "primary"
                  )
                ),
                column(
                  width = 6,
                  box(
                    title = "Expression heatmap by signagure",solidHeader = TRUE,status = "primary"
                  )
                )
              ),
      tabItem("iter")

      # analysis  panel ----
    )
    # Boxes need to be put in a row (or column)


  )
  )
))
