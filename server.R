# server.R
library(data.table)
library(survival)
library(survminer)
library(feather)
library(patchwork)
library(dplyr)

options(shiny.maxRequestSize = 500 * 1024 ^ 2)
#correl.table <- NULL
server <- function(input, output, session) {
  surv_table <- reactiveVal()
  
  
  output$survival_time_label <- renderUI({
    req(input$survFile)
    #tables$surv_table <- read_feather(input$survFile)
    surv_table <<-
      read_feather(input$survFile$datapath) %>% filter(cancer_type_code == "LMS") #BAD
    
    selectInput(
      inputId = "survival_time",
      label = "Select the survival time",
      choices = colnames(surv_table),
      multiple = F
    )
  })
  
  
  output$survival_event_label <- renderUI({
    req(input$survFile)
    selectInput(
      inputId = "survival_event",
      label = "Select the survival event",
      choices = colnames(surv_table),
      multiple = F
    )
  })
  
  output$genes_selection <- renderUI({
    req(input$survFile)
    # TODO: sort genes by frequency of mutation
    selectInput(
      inputId = "genes",
      label = "Select gene(s) to stratify.",
      choices = colnames(surv_table),
      multiple = T
    )
  })
  
  
  # Generate a plot of the data. Also uses the inputs to build the
  # plot label. Note that the dependencies on both the inputs and
  # the 'data' reactive expression are both tracked, and all expressions
  # are called in the sequence implied by the dependency graph
  output$plot1 <- renderPlot({
    req(surv_table, input$survival_time, input$survival_event)
    if (length(input$genes) == 0) {
      strata <- rep(1, times = nrow(surv_table))
    } else {
      strata <- surv_table[[input$genes[[1]]]]
    }
    form <-
      do.call(Surv, list(time = surv_table[[input$survival_time]], event = surv_table[[input$survival_event]]))
    print(form)
    print(strata)
    fit <-
      do.call(survfit, list(formula = form ~ strata, data = surv_table))
    
    p <-
      ggsurvplot(
        data = surv_table,
        fit = fit,
        pval = T,
        title = input$genes[[1]],
        legend.labs = c(paste(input$genes[[1]], "wt"), paste(input$genes[[1]], "mut"))
      )
    p
  })
  
  output$plot2 <- renderPlot({
    strata <- surv_table[[input$genes[[2]]]]
    
    form <-
      do.call(Surv, list(time = surv_table[[input$survival_time]], event = surv_table[[input$survival_event]]))
    print(form)
    print(strata)
    fit <-
      do.call(survfit, list(formula = form ~ strata, data = surv_table))
    
    p <-
      ggsurvplot(
        data = surv_table,
        fit = fit,
        pval = T,
        title = input$genes[[2]],
        legend.labs = c(paste(input$genes[[2]], "wt"), paste(input$genes[[2]], "mut"))
      )
    p
  })
  
  
  output$plot3 <- renderPlot({
    # Both WT
    strata <-
      !surv_table[[input$genes[[1]]]] & !surv_table[[input$genes[[2]]]]
    
    form <-
      do.call(Surv, list(time = surv_table[[input$survival_time]], event = surv_table[[input$survival_event]]))
    print(form)
    print(strata)
    fit <-
      do.call(survfit, list(formula = form ~ strata, data = surv_table))
    
    p <-
      ggsurvplot(
        data = surv_table,
        fit = fit,
        pval = T,
        title = paste(input$genes[[1]], " and ", input$genes[[2]], " wt", sep = ""),
        legend.labs = c("At least one mut", "Both wt"),
        risk.table = F
      )
    p
  })
  
  output$plot4 <- renderPlot({
    # Gene 1 mut Gene 2 wt
    strata <-
      surv_table[[input$genes[[1]]]] & !surv_table[[input$genes[[2]]]]
    
    form <-
      do.call(Surv, list(time = surv_table[[input$survival_time]], event = surv_table[[input$survival_event]]))
    print(form)
    print(strata)
    fit <-
      do.call(survfit, list(formula = form ~ strata, data = surv_table))
    
    p <-
      ggsurvplot(
        data = surv_table,
        fit = fit,
        pval = T,
        title = paste(input$genes[[1]], " mut and ", input$genes[[2]], " wt", sep = ""),
        legend.labs = c(
          paste(input$genes[[1]], "wt OR", input$genes[[2]], "mut"),
          paste(input$genes[[1]], "mut and", input$genes[[2]], "wt")
        ),
        risk.table = F
      )
    p
  })
  

  output$plot5 <- renderPlot({
    # Gene 1 wt Gene 2 mut
    strata <-
      !surv_table[[input$genes[[1]]]] & surv_table[[input$genes[[2]]]]
    
    form <-
      do.call(Surv, list(time = surv_table[[input$survival_time]], event = surv_table[[input$survival_event]]))
    print(form)
    print(strata)
    fit <-
      do.call(survfit, list(formula = form ~ strata, data = surv_table))
    
    p <-
      ggsurvplot(
        data = surv_table,
        fit = fit,
        pval = T,
        title = paste(input$genes[[1]], " wt and ", input$genes[[2]], " mut", sep = ""),
        legend.labs = c(
          paste(input$genes[[1]], "mut OR", input$genes[[2]], "wt"),
          paste(input$genes[[1]], "wt and", input$genes[[2]], "mut")
        ),
        risk.table = F
      )
    p
  })

  output$plot6 <- renderPlot({
    # Gene 1 mut Gene 2 mut
    strata <-
      surv_table[[input$genes[[1]]]] & surv_table[[input$genes[[2]]]]
    
    form <-
      do.call(Surv, list(time = surv_table[[input$survival_time]], event = surv_table[[input$survival_event]]))
    print(form)
    print(strata)
    fit <-
      do.call(survfit, list(formula = form ~ strata, data = surv_table))
    
    p <-
      ggsurvplot(
        data = surv_table,
        fit = fit,
        pval = T,
        title = paste(input$genes[[1]], " and ", input$genes[[2]], " mut", sep = ""),
        legend.labs = c("Either wt", "Both mut"),
        risk.table = F
      )
    p
  })
  
  output$plotOutput <- renderUI({
    if (is.null(input$genes)) {
      plotOutput("plot1")
    } else if (length(input$genes) == 1) {
      plotOutput("plot1")
    } else if (length(input$genes) == 2) {
      fluidPage(fluidRow(column(6, plotOutput("plot1")), column(6, plotOutput("plot2"))),
                fluidRow(
                  column(3, plotOutput("plot3")),
                  column(3, plotOutput("plot4")),
                  column(3, plotOutput("plot5")),
                  column(3, plotOutput("plot6"))
                ))
    }
  })
  
  # Generate a summary of the data
  # Could put linear model info here
  # output$summary <- renderPrint({
  #   summary(data())
  # })
  
}