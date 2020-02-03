library(tidyverse)

LB_data_only <- function(px1z1, px1z0, pz1, py1x0, py0x1){
  ((px1z1 - px1z0)*pz1 + px1z0) * (py1x0 - py0x1) - py1x0
}

LB_A1_A2 <- function(py1x1, py1x0, px1z0, px1z1){
  p1 <- -(1-py1x1)*px1z0 - py1x0*(1-px1z0)
  p2 <- -(1-py1x1)*px1z1 - py1x0*(1-px1z1)
  p3 <- py1x1*(px1z0 - px1z1) - py1x0*(1-px1z1) - (1-py1x1)*px1z1
  p4 <- py1x1*(px1z1 - px1z0) - py1x0*(1-px1z0) - (1-py1x1)*px1z0

  return(c(p1=p1,p2=p2,p3=p3,p4=p4))
}

LB_A1_A2_orig <- function(py0x1z0, py1x0z0,
                          py0x1z1, py1x0z1,
                          py1z0, py1z1){
  p1 <- -py0x1z0 - py1x0z0
  p2 <- -py0x1z0 - py1x0z1
  p3 <- py1z0 - py1z1 - py1x0z0 - py0x1z1
  p4 <- py1z1 - py1z0 - py1x0z1 - py0x1z0

  return(c(p1,p2,p3,p4))
}

LB_A3_A4 <- function(py1x1, py1x0, px1z0, px1z1){

  first_four <- LB_A1_A2(py1x1, py1x0, px1z0, px1z1)

  p5 <- py1x1*(px1z0 - px1z1) - py1x0*((1 - px1z1) + (1 - px1z0)) - (1 - py1x1)*px1z1
  p6 <- py1x1*(px1z1 - px1z0) - py1x0*((1 - px1z0) + (1 - px1z1)) - (1 - py1x1)*px1z0
  p7 <- (1 - py1x0)*(px1z1 - px1z0) - (1 - py1x1)*(px1z0 + px1z1) - py1x0*(1 - px1z0)
  p8 <- (1 - py1x0)*(px1z0 - px1z1) - (1 - py1x1)*(px1z1 + px1z0) - py1x0*(1 - px1z1)

  return(c(first_four, p5, p6, p7, p8))

}


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

UB_tri <- function(py1x1, py1x0, px1z){

  if(is.list(px1z))
    px1z <- unlist(px1z)

  ub <- g(i = "0", j = "0",
        py1x1 = py1x1, py1x0 = py1x0,
        px1z = c(px1z0, px1z1)) +
    g(i = "1", j = "1",
      py1x1 = py1x1, py1x0 = py1x0,
      px1z = px1z) - 1

  return(ub)
}

# g_orig <- function(i, j,
#                    py1x1z0, py1x0z0, py0x1z0,
#                    py1x1z1, py1x0z1, py0x1z1,
#                    px1z0, px1z1){
#
#   py0x0z0 <- 1 - py1x1z0 - py1x0z0 - py0x1z0
#   py0x0z1 <- 1 - py1x1z1 - py1x0z1 - py0x1z1
#
#   px0z0 <- 1 - px1z0
#   px0z1 <- 1 - px1z1
#
#   pyxz <- array(data = c(py0x0z0, py1x0z0,
#                          py0x1z0, py1x1z0,
#                          py0x0z1, py1x0z1,
#                          py0x1z1, py1x1z1),
#                 dim = c(2,2,2),
#                 dimnames = list(y = 0:1,x = 0:1, z = 0:1))
#
#   pxz <- array(data = c(px0z0, px1z0,
#                         px0z1, px1z1),
#                dim = c(2,2),
#                dimnames = list(x = 0:1, z = 0:1))
#
#   expr1 <- pyxz[j,i,] + pxz[as.character(1-as.numeric(i)),]
#
#   expr2 <- c(pyxz[j,i,"0"] + pyxz["0", as.character(1-as.numeric(i)), "0"] +
#                pyxz[j,i,"1"] + pyxz["1", as.character(1-as.numeric(i)), "1"],
#              pyxz[j,i,"1"] + pyxz["0", as.character(1-as.numeric(i)), "1"] +
#                pyxz[j,i,"0"] + pyxz["1", as.character(1-as.numeric(i)), "0"])
#
#   return(min(c(min(expr1), min(expr2))))
# }

if(FALSE){
  1 - g_orig(i = "1", j = "0",
             py1x1z0 = py1x1z0,
             py1x0z0 = py1x0z0,
             py0x1z0 = py0x1z0,
             py1x1z1 = py1x1z1,
             py1x0z1 = py1x0z1,
             py0x1z1 = py0x1z1,
             px1z0 = px1z0,
             px1z1 = px1z1) -
    g_orig(i = "0", j = "1",
           py1x1z0 = py1x1z0,
           py1x0z0 = py1x0z0,
           py0x1z0 = py0x1z0,
           py1x1z1 = py1x1z1,
           py1x0z1 = py1x0z1,
           py0x1z1 = py0x1z1,
           px1z0 = px1z0,
           px1z1 = px1z1)

  1 - g(i = "1", j = "0",
        py1x1 = py1x1, py1x0 = py1x0,
        px1z = c(px1z0, px1z1)) -
    g(i = "0", j = "1",
      py1x1 = py1x1, py1x0 = py1x0,
      px1z = c(px1z0, px1z1))

  max(LB_A1_A2(py1x1 = py1x1, py1x0 = py1x0, px1z0 = px1z0, px1z1 = px1z1))

  g(i = "0", j = "0",
    py1x1 = py1x1, py1x0 = py1x0,
    px1z = c(px1z0, px1z1)) +
    g(i = "1", j = "1",
      py1x1 = py1x1, py1x0 = py1x0,
      px1z = c(px1z0, px1z1)) - 1

  bpbounds::bpbounds(tabp_from_pyx_pxz(pyx, pxz))
}

data_for_plots <- expand_grid(px1z1 = seq(0, 1, by = 0.1),
            px1z0 = seq(0, 1, by = 0.1),
            pz1 = c(0.2, 0.4, 0.8),
            py1x0 = c(0, 0.2, 0.4, 0.6, 0.8, 1),
            py0x1 = c(0, 0.2, 0.4, 0.6, 0.8, 1)) %>%
  mutate(strength = px1z1 - px1z0,
         LB_data_only = LB_data_only(px1z1, px1z0, pz1, py1x0, py0x1),
         LB_A1_A2 = pmap(.l = list(py1x1 = 1-py0x1, py1x0 = py1x0, px1z0 = px1z0, px1z1 = px1z1),
                         .f = LB_A1_A2),
         LB_A3_A4 = pmap(.l = list(py1x1 = 1-py0x1, py1x0 = py1x0, px1z0 = px1z0, px1z1 = px1z1),
                         .f = LB_A3_A4),
         `py1x0 - py0x1` = factor(paste(py1x0, py0x1, sep = ' - ')),
         y_diffs = py1x0 - py0x1)

y_diffs_neg <- data_for_plots %>%
  filter(py0x1 > py1x0) %>%
  #arrange(`py1x0 - py0x1`) %>%
  ggplot(aes(x = strength, y = LB_data_only,
             group = px1z0, col = as.character(px1z0))) +
    geom_line() +
    facet_grid(py0x1 + py1x0 ~ pz1,
               labeller = labeller(.rows = label_both, .cols = label_both),
               scales = "free_y") +
    guides(color = guide_legend(title = "px1z0",
                                nrow = 1,
                                title.position = "top")) +
    theme_bw() +
    theme(legend.position = "top")

y_diffs_pos <- data_for_plots %>%
  filter(py0x1 < py1x0) %>%
  #arrange(`py1x0 - py0x1`) %>%
  ggplot(aes(x = strength, y = LB_data_only,
             group = px1z0, col = as.character(px1z0))) +
  geom_line() +
  facet_grid(py0x1 + py1x0 ~ pz1,
             labeller = labeller(.rows = label_both, .cols = label_both),
             scales = "free_y") +
  guides(color = guide_legend(title = "px1z0",
                              nrow = 1,
                              title.position = "top")) +
  theme_bw() +
  theme(legend.position = "top")

y_diffs_0 <- data_for_plots %>%
  filter(py1x0 == py0x1) %>%
  ggplot(aes(x = strength, y = LB_data_only,
             group = px1z0, col = as.character(px1z0))) +
    geom_line() +
    facet_grid(py1x0 + py0x1 ~ pz1,
               labeller = labeller(.rows = label_both, .cols = label_both),
               scales = "free_y") +
    guides(color = guide_legend(title = "px1z0",
                                nrow = 1,
                                title.position = "top")) +
    theme_bw() +
    theme(legend.position = "top")





# y_diffs
# no_y_diffs
