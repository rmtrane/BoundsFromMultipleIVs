# Nonparametric Bounds on ATE using Simulations

Here we will explore the effect of instrumental variables on nonparameteric bounds on the average treatment effect (ATE) using simulations. The main product is the Shiny application that let's you specify a setup, and with the push of a button, simulate data, and visualize the results.

## How to run the Shiny app

Make sure you have the following packages using the code below:

```
install.packages(c('shiny','bpbounds','glue','dagitty','ggdag','ggrepel','tidyverse'))
```

Once that is done, you can run the application using this line of code:

```
shiny::runGitHub("rmtrane/BoundsFromMultipleIVs", subdir = "shiny_simulation")
```
