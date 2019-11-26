## Run bpbounds on multiple IVs
my_bpbounds <- function(input_data,
                        Y = "Y", X = "X", Z = NULL, ...){

  if(is.null(Z)){
    Z <- str_extract_all(colnames(input_data), pattern = "Z[0-9]+") %>% unlist()
  }

  res <- tibble(IV = Z,
                bpbounds = map(IV, function(Z, ...){

                  tabs_formula <- as.formula(glue::glue("~ {X} + {Y} + {Z}"))
                  tabled_data <- xtabs(tabs_formula, data = input_data)
                  return(bpbounds(tabled_data, ...))
                }))
  return(res)
}

## Tidy up output from bpbounds
tidy_bpbounds <- function(bpbs){
  bpbs_summary <- summary(bpbs)

  if(bpbs_summary$inequality){
    output <- bpbs_summary$bounds %>%
      mutate_if(is.factor, as.character) %>%
      mutate(Assumption = "IV Inequality",
             `Assumption Holds` = bpbs_summary$inequality)
  } else {
    output <- tibble(Assumption = "IV Inequality",
                     `Assumption Holds` = bpbs_summary$inequality,
                     `Causal parameter` = NA_character_,
                     `Lower bound` = NA_real_,
                     `Upper bound` = NA_real_)
  }

  if(bpbs_summary$monoinequality){
    output <- bpbs_summary$bounds %>%
      mutate_if(is.factor, as.character) %>%
      mutate(Assumption = "Monotonicity",
             `Assumption Holds` = bpbs_summary$monoinequality) %>%
      bind_rows(output)
  } else {
    output <- bind_rows(output,
                        tibble(Assumption = "Monotonicity",
                               `Assumption Holds` = bpbs_summary$monoinequality,
                              `Causal parameter` = NA_character_,
                              `Lower bound` = NA_real_,
                              `Upper bound` = NA_real_))
  }

  out <- output %>%
    dplyr::select(Assumption, `Assumption Holds`, everything()) %>%
    arrange(Assumption)

  return(as_tibble(out))
}


## Simulate from bivariate normal
rbivNorm <- function(n, mu = 0, var1 = 1, var2 = 1, rho = 0, probs = c(0.25, 0.5, 0.25)) {
  tmp <- MASS::mvrnorm(n = n,
                       mu = rep(mu, 2),
                       Sigma = matrix(data = c(var1, rho*sqrt(var1)*sqrt(var2),
                                               rho*sqrt(var1)*sqrt(var2), var2),
                                      nrow = 2, byrow = FALSE))

  tmp <- matrix(tmp, nrow = n)
  out <- apply(tmp, 2,
               function(x){
                 as.numeric(cut(x, breaks = qnorm(unique(c(0, cumsum(probs), 1))))) - 1
               })

  return(out)
}

## Simulate data manually (simcausal not working at all...)
simulate_data <- function(sample_size = 10000,
                          indIVs_on_X = NULL, IVs_ps = c(0.5, 0.5), indIVs_on_Y = NULL,
                          X_intercept = 0.05, U_on_X = 0.1,
                          Y_intercept = 0.05, X_on_Y = 0.05, U_on_Y = 0.1,
                          rho = NULL, depIVs_on_X = NULL, dep_probs = c(0.5, 0.5), depIVs_on_Y = NULL){

  ## Default IV effects on Y to 0
  if(is.null(indIVs_on_Y))
    indIVs_on_Y <- rep(0, length(indIVs_on_X))

  if(is.null(depIVs_on_Y))
    depIVs_on_Y <- rep(0, length(depIVs_on_X))

  ## Start tibble with just why. Include effects on X and Y for later use.
  final_data <- tibble(var = "U",
                       var_on_X = U_on_X,
                       var_on_Y = U_on_Y,
                       simulated = map(var, function(x){
                         as.numeric(rbernoulli(n = sample_size, p = 0.5))
                       }))

  ## If independent IVs are included...
  if(length(indIVs_on_X) > 0){
    ## ... create tible with simulations of these. Include effects on X and Y for later use.
    out_indIVs <- tibble(var = paste0("Z", seq_along(indIVs_on_X)),
                         var_on_X = indIVs_on_X,
                         var_on_Y = indIVs_on_Y,
                         simulated = map(var_on_X, function(x){
                           ## Simulated IVs by sampling from 1,2 (or 1,2,3) with replacement.
                           ## Prob of each is given as argument IVs_ps. Subtract one to turn into
                           ## 0,1 (or 0,1,2) variable.
                           sample(1:length(IVs_ps), replace = TRUE, size = sample_size,
                                  prob = IVs_ps) - 1
                         }))

    ## Add to full data
    final_data <- bind_rows(
      final_data,
      out_indIVs
    )

  }

  ## If idependent IVs are included...
  if(length(depIVs_on_X) > 0){
    tmp <- rbivNorm(n = sample_size,
                    rho = rho, probs = dep_probs)

    colnames(tmp) <- paste0("Z", length(indIVs_on_X) + c(1,2))
    out_depIVs <- as_tibble(tmp) %>%
      gather(key = var, value = simulated) %>%
      nest(simulated = simulated) %>%
      mutate(var_on_X = depIVs_on_X,
             var_on_Y = depIVs_on_Y,
             simulated = map(simulated, function(x) pull(x, simulated)))

    final_data <- bind_rows(
      final_data,
      out_depIVs
    )
  }

  ## Some manipulations
  tmp_data <- final_data %>%
    ## Multiple vars by effects on X and Y, respectively
    mutate(weighted_for_X = map2(var_on_X, simulated,
                                 function(x,y){
                                   tibble(weights_for_X = x*y)
                                 }),
           weighted_for_Y = map2(var_on_Y, simulated,
                                 function(x,y){
                                   tibble(weights_for_Y = x*y)
                                 })) %>%
    ## Expand simulated and weighted observations
    unnest(cols = c(simulated, weighted_for_X, weighted_for_Y)) %>%
    ## Group by var, and create id column
    group_by(var) %>%
    mutate(i = row_number()) %>%
    ## Group by id, and get contribution from U and Z's on X and Y
    group_by(i) %>%
    mutate(pX = sum(weights_for_X) + X_intercept,
           pY = sum(weights_for_Y) + Y_intercept) %>%
    ungroup() %>%
    select(var, simulated, pX, pY, i)

  ## Simulate X
  out_data <- tmp_data %>%
    pivot_wider(id_cols = c(i, var, pX, pY), names_from = var, values_from = simulated) %>%
    mutate(X = map_dbl(pX, ~rbinom(n = 1, size = 1, prob = min(.x, 1))),
           pY = pY + X*X_on_Y,
           Y = map_dbl(pY, ~rbinom(n = 1, size = 1, prob = min(.x, 1)))) %>%
    select(-i)

  #message(paste("Number of pX > 1:", sum(out_data$pX > 1)))
  #message(paste("Number of pY > 1:", sum(out_data$pY > 1)))

  coefs <- tibble(
    coef = c(X_intercept,
             indIVs_on_X,
             depIVs_on_X,
             U_on_X,
             Y_intercept,
             X_on_Y,
             U_on_Y),
    effect = c("Xintercept",
               if(length(indIVs_on_X) > 0){
                 glue("Z{1:length(indIVs_on_X)}_on_X")
               } else {
                 NULL
               },
               if(length(depIVs_on_X) > 0){
                 glue("Z{length(indIVs_on_X) + c(1,2)}_on_X")
               } else {
                 NULL
               },
               "U_on_X",
               "Yintercept",
               "X_on_Y",
               "U_on_Y")
  )

  if(length(depIVs_on_X) > 0){
    coefs <- coefs %>%
      bind_rows(
        tibble(coef = rho,
               effect = glue("Z{length(indIVs_on_X)+1}_on_Z{length(indIVs_on_X)+2}"))
      )
  }

  if(any(indIVs_on_Y > 0)){
    coefs <- bind_rows(
      coefs,
      tibble(coef = indIVs_on_Y,
             effect = glue("Z{1:length(indIVs_on_Y)}_on_Y"))
    ) %>%
      filter(coef != 0)
  }

  if(any(depIVs_on_Y > 0)){
    coefs <- bind_rows(
      coefs,
      tibble(coef = depIVs_on_Y,
             effect = glue("Z{length(depIVs_on_Y) + c(1,2)}_on_Y"))
    ) %>%
      filter(coef != 0)
  }

  return(list(simulated_data = select(out_data, -pX, -pY),
              coefficients = coefs,
              IVs_ps = IVs_ps,
              dep_probs = dep_probs))
}

from_simulated_to_dag <- function(simulated_data, layout = "auto"){

  coefs <- simulated_data$coefficients

  tmp <- coefs %>%
    filter(str_detect(effect, "_on_")) %>%
    mutate(edges = str_replace_all(effect, pattern = "_on_", "->") %>% paste0(., ";"))

  all_nodes <- tmp$effect %>% str_split("_on_", simplify = T) %>% .[,1] %>% unique()

  DAG <- paste0("dag{",
                paste(tmp$edges, collapse = ""),
                "}") %>%
    dagitty()

  tidy_dag <- tidy_dagitty(DAG,
                           layout = layout)

  tidy_dag$data <- tidy_dag$data %>%
    left_join(
      coefs %>%
        separate(effect, into = c("name", "to"), sep = "_on_", fill = "right") %>%
        filter(!is.na(to))
    )

  any_corrZs <- str_detect(coefs$effect,
                           pattern = "Z[0-9]_on_Z[0-9]")

  if(any(any_corrZs)){
    corrZs <- str_split(coefs$effect[any_corrZs], pattern = "_on_")[[1]]

    tidy_dag$data <- tidy_dag$data %>%
      mutate(direction = if_else(name %in% corrZs & to %in% corrZs,
                                 factor("<->", levels = levels(direction)),
                                 direction))
  }

  tidy_dag$data <- tidy_dag$data %>%
    mutate(label_x = (x + xend)/2,
           label_y = (y + yend)/2)

  return(tidy_dag)
}


plot_bounds <- function(my_bds, sim_data_coefficients){

  tidied_output <- my_bds %>%
    mutate(tidy_output = map(bpbounds, tidy_bpbounds),
           `n Z categories` = map_dbl(bpbounds, function(x) x$nzcats)) %>%
    unnest(tidy_output)

  for_plot <- tidied_output %>%
    filter(`Assumption Holds`) %>%
    filter(Assumption == "IV Inequality", `Causal parameter` == "ACE") %>%
    select(-bpbounds, -Assumption, -`Assumption Holds`, -`Causal parameter`) %>%
    left_join(
      sim_data_coefficients %>%
        filter(str_detect(effect, "Z[0-9]+_on_X")) %>%
        separate(effect, into = c("IV", "X"), sep = "_on_") %>%
        select(-X)
    )

  IVs <- for_plot$IV
  IV_coefs <- for_plot$coef

  true_effect <- sim_data_coefficients %>% filter(effect == "X_on_Y") %>% pull(coef)

  ggplot(for_plot,
         aes(ymin = `Lower bound`, ymax = `Upper bound`, x = seq_along(IV))) +
    geom_hline(yintercept = true_effect, color = "red", linetype = "dashed") +
    geom_errorbar(width = 0.25) +
    scale_x_continuous(name = "IV",
                       labels = IVs,
                       breaks = seq_along(IVs),
                       sec.axis = sec_axis(trans = function(x) x,
                                           breaks = seq_along(IVs),
                                           labels = IV_coefs,
                                           name = "Effect on X")) +
    scale_y_continuous(limits = c(-1,1)) +
    theme_bw()
}