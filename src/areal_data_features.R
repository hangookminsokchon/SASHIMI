######################################
###### Areal Feature Extraction ######
######################################
# Input: 
#   T, S, L - ppp objects (Tumor, Stromal, Lymphocyte)
# Output: 
#   Data frames with scalar areal features per cell type

###### Moran's I ######
moran_I <- function(T, S, L) {
  nb <- cell2nb(20, 20, type = "queen")
  weight <- nb2mat(nb, style = "B", zero.policy = TRUE)
  
  mi_t <- Moran.I(as.vector(T), weight)
  mi_s <- Moran.I(as.vector(S), weight)
  mi_l <- Moran.I(as.vector(L), weight)
  
  data.frame("MI.T" = mi_t, "MI.S" = mi_s, "MI.I" = mi_l)
}

###### Morisita-Horn Index ######
morisita_horn <- function(T, S, L) {
  compute_mh <- function(x, y) {
    p1 <- x / sum(x)
    p2 <- y / sum(y)
    (2 * sum(p1 * p2)) / (sum(p1^2) + sum(p2^2))
  }
  
  mh_ts <- compute_mh(as.vector(T), as.vector(S))
  mh_ti <- compute_mh(as.vector(T), as.vector(L))
  mh_is <- compute_mh(as.vector(L), as.vector(S))
  
  data.frame("MH.TS" = mh_ts, "MH.TI" = mh_ti, "MH.IS" = mh_is)
}

###### Lee's L ######
lee_L <- function(T, S, L) {
  nb <- cell2nb(20, 20, type = "queen")
  lw <- nb2listw(nb, style = "W")
  
  l_ti <- lee.test(as.vector(T), as.vector(L), listw = lw)$estimate[1]
  l_ts <- lee.test(as.vector(T), as.vector(S), listw = lw)$estimate[1]
  l_is <- lee.test(as.vector(L), as.vector(S), listw = lw)$estimate[1]
  
  data.frame("Lee.TI" = l_ti, "Lee.TS" = l_ts, "Lee.IS" = l_is)
}

###### Geary's C ######
geary_C <- function(T, S, L) {
  nb <- cell2nb(20, 20, type = "queen")
  lw <- nb2listw(nb, style = "C")
  
  gc_t <- as.numeric(geary.test(T, listw = lw)$estimate[1])
  gc_s <- as.numeric(geary.test(S, listw = lw)$estimate[1])
  gc_l <- as.numeric(geary.test(L, listw = lw)$estimate[1])
  
  data.frame("GearyC.T" = gc_t, "GearyC.S" = gc_s, "GearyC.I" = gc_l)
}

###### Bhattacharyya Coefficient ######
bhattacharyya <- function(T, S, L) {
  bc_ts <- bhatt.coeff(as.numeric(T), as.numeric(S))
  bc_ti <- bhatt.coeff(as.numeric(T), as.numeric(L))
  bc_is <- bhatt.coeff(as.numeric(L), as.numeric(S))
  
  data.frame("BC.TS" = bc_ts, "BC.TI" = bc_ti, "BC.IS" = bc_is)
}

###### Placeholders for future features ######
# clark_evans <- function(...) { ... }
# quadrat_statistic <- function(...) { ... }
# join_count <- function(...) { ... }
# jaccard_index <- function(...) { ... }
# dice_sorensen_index <- function(...) { ... }










