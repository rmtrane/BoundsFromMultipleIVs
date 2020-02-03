library(tidyverse)

g <- function(i, j, py1x1, py1x0, px1z){
  # i = {0,1}, j = {0,1}}
  # py1x0: P(Y = 1 | X = 0)
  # py1x1: P(Y = 1 | X = 1)
  # px1z: 2 x k matrix with P(X = i | Z = j) as entry (i,j)

  pyx <- tribble(
    ~y, ~x, ~p,
    0, 0, 1-py1x0,
    0, 1, 1-py1x1,
    1, 0, py1x0,
    1, 1, py1x1
  ) %>%
    pivot_wider(names_from = x, values_from = p) %>%
    column_to_rownames("y") %>%
    as.matrix()



  pxz <- expand_grid(x = c(0,1),
                     z = 1:length(px1z) - 1) %>%
    mutate(p = case_when(x == 0 ~ 1-px1z[z+1],
                         x == 1 ~ px1z[z+1],
                         TRUE ~ NA_real_)) %>%
    pivot_wider(names_from = z, values_from = p) %>%
    column_to_rownames("x") %>%
    as.matrix()

  expr1 <- min(pyx[j, i]*pxz[i, ] + pxz[as.character(1 - as.numeric(i)), ])

  # expr2 <- pyx[1, 1 - i + 1]*outer(pxz[1 - i + 1, ], pxz[1 - i + 1, ], `-`) +
  #   pyx[j + 1, i + 1]*outer(pxz[i + 1,], pxz[i + 1,], `+`) +
  #   matrix(rep(pxz[1 - i + 1,], ncol(pxz)),
  #          nrow = ncol(pxz), byrow = TRUE)

  expr2 <- min(c(pyx["1",as.character(1-as.numeric(i))]*(pxz[as.character(1-as.numeric(i)), "0"] - pxz[as.character(1-as.numeric(i)), "1"]) + pyx[j, i]*(pxz[i, "0"] + pxz[i, "1"]) + pxz[as.character(1-as.numeric(i)), "1"],
                 pyx["1",as.character(1-as.numeric(i))]*(pxz[as.character(1-as.numeric(i)), "1"] - pxz[as.character(1-as.numeric(i)), "0"]) + pyx[j, i]*(pxz[i, "0"] + pxz[i, "1"]) + pxz[as.character(1-as.numeric(i)), "0"]))

  return(min(expr1, expr2))
}

ATE_bounds_tri <- function(py1x1, py1x0, px1z){

  if(is.list(px1z))
    px1z <- unlist(px1z)

  lb <- 1 - g(i = "1", j = "0",
              py1x1 = py1x1, py1x0 = py1x0,
              px1z = px1z) -
    g(i = "0", j = "1",
      py1x1 = py1x1, py1x0 = py1x0,
      px1z = px1z)

  ub <- g(i = "0", j = "0",
          py1x1 = py1x1, py1x0 = py1x0,
          px1z = px1z) +
    g(i = "1", j = "1",
      py1x1 = py1x1, py1x0 = py1x0,
      px1z = px1z) - 1

  return(tibble(LB = lb, UB = ub))
}

library(furrr)
plan(multiprocess)

new_dat <- expand_grid(px1z0 = seq(0, 1, by = 0.1),
                       px1z1 = seq(0, 1, by = 0.1),
                       px1z2 = seq(0, 1, by = 0.1),
                       py1x0 = seq(0, 1, by = 0.2),
                       py1x1 = seq(0, 1, by = 0.2)) %>%
  head(n = 100) %>%
  mutate(bounds = future_pmap(.l = list(a = py1x1, b = py1x0, c = px1z0, d = px1z1, e = px1z2),
                              .f = function(a,b,c,d,e) ATE_bounds_tri(py1x1 = a, py1x0 = b, px1z = c(c,d,e)),
                              .progress = TRUE))

write_rds(x = new_dat, path = "new_dat.Rds")
