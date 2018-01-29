ui <- fluidPage(
  # Application title
  titlePanel("Correlation between Copy Number and Gene Expression"),
  
  # Sidebar with controls to select the random distribution type
  # and number of observations to generate. Note the use of the br()
  # element to introduce extra vertical spacing
  sidebarPanel(width = 2,
    tags$h3("Needed for plot:"),
    fileInput('survFile', 'Upload dataset'),
    uiOutput("survival_time_label"),
    uiOutput("survival_event_label"),
    tags$hr(),
    uiOutput("genes_selection")
    # uiOutput("band"),
    # uiOutput("category")
    
  ),
  
  # Show a tabset that includes a plot, summary, and table view
  # of the generated distribution
  mainPanel(width = 10,
    tabsetPanel(
      tabPanel("Plot",uiOutput("plotOutput"))
      #tabPanel("Summary", verbatimTextOutput("summary")), 
      #tabPanel("Table", dataTableOutput("table"))
    )
  ),
  mainPanel(
    tags$a(href="http://nickgros.com/shiny.html", "nickgros.com")
  )
)