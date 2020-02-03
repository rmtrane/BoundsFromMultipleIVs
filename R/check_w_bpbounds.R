cp <- c(.0064, 0, .9936, 0, .0028, .001, .1972, .799)

cp <- runif(n = 8, 0, 1)

cp[1:4] <- cp[1:4]/sum(cp[1:4])
cp[5:8] <- cp[5:8]/sum(cp[5:8])

# names(cp) <- c("py0x0|z0", "py0x1|z0", "py1x0|z0", "py1x1|z0",
#                "py0x0|z1", "py0x1|z1", "py1x0|z1", "py1x1|z1")

tabp <- as.table(array(
  cp,
  dim = c(2, 2, 2),
  dimnames = list(
    x = c(0, 1),
    y = c(0, 1),
    z = c(0, 1)
  )
))

tab1dat <- data.frame(
  z = c(0, 0, 1, 1, 1, 1, 0, 0),
  x = c(0, 0, 0, 0, 1, 1, 1, 1),
  y = c(0, 1, 0, 1, 0, 1, 0, 1),
  freq = c(74, 11514, 34, 2385, 12, 9663, 0, 0)
)

tab1inddat = uncount(tab1dat, freq)
xt = xtabs(~ x + y + z, data = tab1inddat)
tabp = prop.table(xt, margin = 3)

pyx <- function(i, j, tabp, pz = c(0.5, 0.5)){

  nom <- sum(tabp[i,j,]*pz)
  den <- sum(tabp[i,,]*pz)

  return(nom/den)
}

# py1x1 <- pyx("1","1", tabp = cp_tab$cp_tab, pz = cp_tab$pz$percent)
# py1x0 <- pyx("0","1", tabp = cp_tab$cp_tab, pz = cp_tab$pz$percent)
#
# px1z0 <- sum(cp_tab$cp_tab["1",,"0"])
# px1z1 <- sum(cp_tab$cp_tab["1",,"1"])

py1x1 <- pyx("1","1", tabp = tabp)
py1x0 <- pyx("0","1", tabp = tabp)

px1z0 <- sum(tabp["1",,"0"])
px1z1 <- sum(tabp["1",,"1"])


bpbounds::bpbounds(tabp)

LB_data_only(px1z1 = px1z1,
             px1z0 = px1z0,
             pz1 = 0.5,
             py1x0 = py1x0,
             py0x1 = 1-py1x1)


py0x0z0 = tabp["0", "0", "0"]
py0x0z1 = tabp["0", "0", "1"]
py0x1z0 = tabp["1", "0", "0"]
py0x1z1 = tabp["1", "0", "1"]
py1x0z0 = tabp["0", "1", "0"]
py1x0z1 = tabp["0", "1", "1"]
py1x1z0 = tabp["1", "1", "0"]
py1x1z1 = tabp["1", "1", "1"]
py1z0 = sum(tabp[,"1","0"])
py1z1 = sum(tabp[,"1","1"])

px1z0 = sum(tabp["1",,"0"])
px1z1 = sum(tabp["1",,"1"])

py1x1 = (py1x1z0 + py1x1z1)/(py0x1z0 + py0x1z1 + py1x1z0 + py1x1z1)
py1x0 = (py1x0z0 + py1x0z1)/(py0x0z0 + py0x0z1 + py1x0z0 + py1x0z1)
px1z0 = px1z0
px1z1 = px1z1


LB_A1_A2(py1x1 = (py1x1z0 + py1x1z1)/(py0x1z0 + py0x1z1 + py1x1z0 + py1x1z1),
         py1x0 = (py1x0z0 + py1x0z1)/(py0x0z0 + py0x0z1 + py1x0z0 + py1x0z1),
         px1z0 = px1z0,
         px1z1 = px1z1)



LB_A1_A2_orig(py0x1z0 = tabp["1", "0", "0"],
              py1x0z0 = tabp["0", "1", "0"],
              py0x1z1 = tabp["1", "0", "1"],
              py1x0z1 = tabp["0", "1", "1"],
              py1z0 = sum(tabp[,"1","0"]),
              py1z1 = sum(tabp[,"1","1"]))


LB_A3_A4(py1x1 = py1x1,
         py1x0 = py1x0,
         px1z0 = px1z0,
         px1z1 = px1z1) %>% max()


cps <- function(data, x, y, z){

  x <- enquo(x)
  y <- enquo(y)
  z <- enquo(z)

  my_prop_table <- for_analysis %>%
    tabyl(!!x, !!y, !!z) %>%
    adorn_percentages(denominator = 'all')

  pzs <- for_analysis %>%
    tabyl(!!z)

  tmp_array <- array(dim = c(2,2,2))
  tmp_array[,,1] <- as.matrix(my_prop_table[[1]][,-1])
  tmp_array[,,2] <- as.matrix(my_prop_table[[2]][,-1])

  # dimnames(tmp_array) <- list(my_prop_table[[1]][[1]],
  #                             colnames(my_prop_table[[1]])[-1],
  #                             names(my_prop_table))

  dimnames(tmp_array) <- list(c(0,1), c(0,1), c(0,1))

  return(list(cp_tab = as.table(tmp_array),
              pz = pzs))
}


cp_tab <- cps(data = for_analysis, x = medicaid, y = any_visit_ed, z = lottery)


bpbounds::bpbounds(p = cp_tab$cp_tab)



tabp_from_pyx_pxz <- function(pyx, pxz){

  if(identical(dim(pyx), c(2,2)))
    stop("x and y must both be binary")

  tabp <- array(data = NA,
                dim = c(2, 2, ncol(pxz)),
                dimnames = list(x = 0:1, y = 0:1, z = 0:(ncol(pxz) - 1)))

  for (z in 1:ncol(pxz)){
    tabp[,,z] <- rbind(pyx[,1]*pxz[1,z],
                       pyx[,2]*pxz[2,z])
  }

  return(as.table(tabp))
}





