ui <- navbarPage(
  "Non-parametric bounds on ATE using multiple IVs",
  tabPanel(
    "Simulation",
    column(
      width = 4,
      h2("Simulation Parameters"),
      wellPanel(
        numericInput(inputId = "sample_size",
                     label = "Sample Size",
                     value = 10000,
                     min = 10,
                     max = 100000),
        checkboxInput(inputId = "invalid_IVs",
                      label = "Allow Invalid IVs",
                      value = FALSE),
        checkboxInput(inputId = "include_dep_IVs",
                      label = "Include a pair of dependent variables",
                      value = FALSE)
      ),
      wellPanel(
        h3("Unmeasured Confounder"),
        numericInput(inputId = "pU",
                     label = HTML("P(U = 1)"),
                     value = 0.5),
        column(6, numericInput(inputId = "U_on_X",
                               label = HTML("Effect of U on X (&alpha;<sub>U</sub>)"),
                               value = 1)),
        column(6, numericInput(inputId = "U_on_Y",
                               label = HTML("Effect of U on Y (&beta;<sub>U</sub>)"),
                               value = 1))
      ),
      wellPanel(
        h3("Effect of X"),
        numericInput(inputId = "X_on_Y",
                     label = HTML("Effect of X on Y (&beta;<sub>X</sub>)"),
                     value = 2)
      ),
      wellPanel(
        h3("Instrumental Variables"),
        sliderInput(inputId = "n_indIVs", label = "Number of independent IVs to include",
                    value = 1, min = 0, max = 10, step = 1),
        radioButtons(inputId = "n_cats",
                     label = "Number of categories of Z",
                     choices = c("2", "3"), inline = TRUE),
        uiOutput(outputId = "IVs_ps"),
        h4("Effect of independent IVs on X"),
        uiOutput(outputId = "indIVs_on_X"),
        uiOutput(outputId = "indIVs_on_Y"),
        uiOutput(outputId = "depIVs")
      ),
      actionButton(inputId = "simulate_data", label = "Simulate data")
    ),
    column(
      h2("Simulation Results"),
      conditionalPanel("input.simulate_data == 0",
                       h3("Choose parameters on the left, and click 'Simulate Data' at the bottom of the sidebar panel.")),
      conditionalPanel("input.simulate_data > 0",
                       plotOutput("DAG_plotted")),
      conditionalPanel("input.simulate_data > 0",
                       plotOutput("plot_results")),
      conditionalPanel("input.simulate_data > 0",
                       DT::dataTableOutput("tidy_results")),
      conditionalPanel("input.simulate_data > 0",
                       DT::dataTableOutput("simulated_data")),
      width = 8

    )
  ),
  tabPanel(
    "Documentation",
    #h2("t")
    includeMarkdown("documentation.md"),
    withMathJax()
    #includeHTML("documentation.html")
  )
)
