library(shiny)
library(bpbounds)
library(glue)
library(dagitty)
library(ggdag)
library(ggrepel)
library(tidyverse)

source("../R/misc.R")

server <- function(input, output){
  output$indIVs_on_X <- renderUI({
    if(input$n_indIVs > 0){
      lapply(1:input$n_indIVs,
             function(i){
               numericInput(inputId = paste0("Z", i, "_on_X"),
                            label = HTML(glue("Z<sub>{i}</sub> (&alpha;<sub>i</sub>)")),
                            value = 0.2*i)
             })
    }
  })

  output$indIVs_on_Y <- renderUI({
    if(input$n_indIVs > 0 & input$invalid_IVs){
      lapply(0:input$n_indIVs,
             function(i){
               if(i < 1){
                 h3("Effect of independent IVs on Y")
               } else {
                 numericInput(inputId = paste0("Z", i, "_on_Y"),
                              label = HTML(glue("Z<sub>{i}</sub> (&beta;<sub>{i}</sub>)")),
                              value = 0)
               }
             })
    }
  })

  output$IVs_ps <- renderUI({
    lapply(1:(as.numeric(input$n_cats)-1),
           function(i){
             numericInput(inputId = paste0("IVs_ps", i),
                          label = HTML(glue("P(Z<sub>j</sub> = {i-1}) for all j")),
                          value = round(1/as.numeric(input$n_cats), digits = 3),
                          min = 0,
                          max = 1)
           })
  })

  output$depIVs <- renderUI({
    if(input$include_dep_IVs){
      outs <- list(h3("Correlated Instrumental Variables"),
                   sliderInput(inputId = "rho",
                               label = "Correlation between dependend IVs (&rho;)",
                               min = 0, max = 1, value = 0.5, step = 0.1),
                   h3("Effects of Correlated IVs on X"),
                   numericInput(inputId = "depIV1_on_X",
                                label = HTML("Z<sub>D1</sub> (&alpha;<sub>D1</sub>)"),
                                value = 0.1),
                   numericInput(inputId = "depIV2_on_X",
                                label = HTML("Z<sub>D1</sub> (&alpha;<sub>D2</sub>)"),
                                value = 0.1)
      )

      return(outs)
    }

  })

  RVs <- reactiveValues()

  observeEvent(input$simulate_data, {
    if(input$n_indIVs > 0){
      RVs$indIVs_on_X <- sapply(1:input$n_indIVs,
                                function(i) input[[paste0("Z", i, "_on_X")]])
      RVs$indIVs_on_Y <- sapply(1:input$n_indIVs,
                                function(i) input[[paste0("Z", i, "_on_Y")]])
    } else {
      RVs$indIVs_on_X <- NULL
      RVs$indIVs_on_Y <- NULL
    }

    if(!input$invalid_IVs)
      RVs$indIVs_on_Y <- NULL

    RVs$IVs_ps <- sapply(1:(as.numeric(input$n_cats)-1), function(i) input[[paste0("IVs_ps", i)]])

    print(paste("Before:", paste(RVs$IVs_ps, collapse = ", ")))

    RVs$IVs_ps[length(RVs$IVs_ps) + 1] <- 1 - sum(RVs$IVs_ps)

    print(paste("After:", paste(RVs$IVs_ps, collapse = ", ")))

    if(input$include_dep_IVs){
      RVs$depIVs_on_X <- c(input$depIV1_on_X, input$depIV2_on_X)
    } else {
      RVs$depIVs_on_X <- NULL
    }



    RVs$sim_data <- simulate_data(sample_size = input$sample_size,
                                  indIVs_on_X = RVs$indIVs_on_X,
                                  indIVs_on_Y = RVs$indIVs_on_Y,
                                  IVs_ps = RVs$IVs_ps,
                                  U_on_X = input$U_on_X,
                                  U_on_Y = input$U_on_Y,
                                  X_on_Y = input$X_on_Y,
                                  rho = input$rho,
                                  depIVs_on_X = RVs$depIVs_on_X,
                                  dep_probs = RVs$IVs_ps,
                                  X_intercept = 0,
                                  Y_intercept = 0,
                                  pU = input$pU
    )

    readr::write_rds(RVs$sim_data, path = "sim_data_from_shiny.Rds")
  })

  observeEvent(RVs$sim_data, {
    RVs$multiple_IVs <- my_bpbounds(input_data = RVs$sim_data$simulated_data,
                                    Y = "Y", X = "X",
                                    Z = colnames(RVs$sim_data$simulated_data)[str_detect(colnames(RVs$sim_data$simulated_data), pattern = "(Z.)")])

    output$tidy_results <- DT::renderDT({
      RVs$multiple_IVs %>%
        mutate(tidy_output = map(bpbounds, tidy_bpbounds),
               `n Z categories` = map_dbl(bpbounds, function(x) x$nzcats)) %>%
        unnest(tidy_output) %>%
        DT::datatable(filter = "top")
    })

    output$simulated_data <- DT::renderDT({
      DT::datatable(RVs$sim_data$simulated_data,
                    filter = "top")
    })

    output$DAG_plotted <- renderPlot({
      tidyDAG <- from_simulated_to_dag(RVs$sim_data)

      ggdag(tidyDAG) +
        geom_label_repel(aes(label = coef,
                             x = label_x,
                             y = label_y)) +
        theme_dag_blank() +
        ggtitle("DAG structure used for simulation")
    })

    output$plot_results <- renderPlot({
      plot_bounds(RVs$multiple_IVs, RVs$sim_data$coefficients)
    })

  })

}
