# Description: Compute a suite of scalar spatial statistics (grid-based) 
#   for Tumor (T), Stromal (S), and Lymphocyte/Immune (L) point patterns.
#
# Dependencies: spatstat, spdep, dispRity, dplyr
#
# All functions use a fixed 20×20 grid and Queen‐type adjacency.
#-------------------------------------------------------------------------------

library(spatstat)        # ppp, quadratcount, superimpose, clarkevans
library(spdep)           # cell2nb, nb2mat, nb2listw, Moran.I, lee.test, geary.test
library(dispRity)        # bhatt.coeff
library(dplyr)           # tapply, etc.

#-------------------------------------------------------------------------------
#' Moran's I for Tumor, Stroma, and Immune
#'
#' Compute Moran's I observed statistic on quadrat counts for three point patterns.
#'
#' @param T A `ppp` object of tumor points.
#' @param S A `ppp` object of stromal points.
#' @param L A `ppp` object of lymphocyte/immune points.
#' @return A 1-row `data.frame` with columns `MI.T`, `MI.S`, `MI.I`.
#' @examples
#' moran_I(tumor_ppp, stroma_ppp, lymph_ppp)
moran_I <- function(T, S, L) {
  # build queen‐adjacency weight matrix for 20×20 grid
  nb     <- cell2nb(20, 20, type = "queen")
  weight <- nb2mat(nb, style = "B", zero.policy = TRUE)
  
  # extract observed Moran's I
  mi_t <- Moran.I(as.vector(quadratcount(T, 20, 20)), weight)$observed
  mi_s <- Moran.I(as.vector(quadratcount(S, 20, 20)), weight)$observed
  mi_l <- Moran.I(as.vector(quadratcount(L, 20, 20)), weight)$observed
  
  data.frame(MI.T = mi_t, MI.S = mi_s, MI.I = mi_l)
}

#-------------------------------------------------------------------------------
#' Morisita–Horn Index for Pairwise Patterns
#'
#' Compute the Morisita–Horn similarity index between each pair of quadrat counts.
#'
#' @param T A `ppp` object of tumor points.
#' @param S A `ppp` object of stromal points.
#' @param L A `ppp` object of lymphocyte/immune points.
#' @return A 1-row `data.frame` with `MH.TS`, `MH.TI`, `MH.IS`.
#' @examples
#' morisita_horn(tumor_ppp, stroma_ppp, lymph_ppp)
morisita_horn <- function(T, S, L) {
  # helper to get flattened counts vector
  get_counts <- function(x) {
    m <- if (inherits(x, "ppp")) quadratcount(x, 20, 20) else x
    as.vector(as.table(m))
  }
  compute_mh <- function(x, y) {
    p1 <- x / sum(x); p2 <- y / sum(y)
    2 * sum(p1 * p2) / (sum(p1^2) + sum(p2^2))
  }
  
  cts_T <- get_counts(T); cts_S <- get_counts(S); cts_L <- get_counts(L)
  mh_ts <- compute_mh(cts_T, cts_S)
  mh_ti <- compute_mh(cts_T, cts_L)
  mh_is <- compute_mh(cts_L, cts_S)
  
  data.frame(MH.TS = mh_ts, MH.TI = mh_ti, MH.IS = mh_is)
}

#-------------------------------------------------------------------------------
#' Lee's L Cross‐Type Autocorrelation
#'
#' Compute Lee’s L statistic for each pair of patterns.
#'
#' @param T A `ppp` object of tumor points.
#' @param S A `ppp` object of stromal points.
#' @param L A `ppp` object of lymphocyte/immune points.
#' @return A 1-row `data.frame` with `Lee.TI`, `Lee.TS`, `Lee.IS`.
#' @examples
#' lee_L(tumor_ppp, stroma_ppp, lymph_ppp)
lee_L <- function(T, S, L) {
  nb  <- cell2nb(20, 20, type = "queen")
  lw  <- nb2listw(nb, style = "W", zero.policy = TRUE)
  
  # flatten counts then apply lee.test
  vT <- as.vector(quadratcount(T, 20, 20))
  vS <- as.vector(quadratcount(S, 20, 20))
  vL <- as.vector(quadratcount(L, 20, 20))
  
  l_ti <- lee.test(vT, vL, listw = lw)$estimate[1]
  l_ts <- lee.test(vT, vS, listw = lw)$estimate[1]
  l_is <- lee.test(vL, vS, listw = lw)$estimate[1]
  
  data.frame(Lee.TI = l_ti, Lee.TS = l_ts, Lee.IS = l_is)
}

#-------------------------------------------------------------------------------
#' Geary's C for Each Pattern
#'
#' Compute Geary's C statistic on quadrat counts.
#'
#' @param T A `ppp` object of tumor points.
#' @param S A `ppp` object of stromal points.
#' @param L A `ppp` object of lymphocyte/immune points.
#' @return A 1-row `data.frame` with `GearyC.T`, `GearyC.S`, `GearyC.I`.
#' @examples
#' geary_C(tumor_ppp, stroma_ppp, lymph_ppp)
geary_C <- function(T, S, L) {
  nb <- cell2nb(20, 20, type = "queen")
  lw <- nb2listw(nb, style = "C", zero.policy = TRUE)
  
  vT <- as.vector(quadratcount(T, 20, 20))
  vS <- as.vector(quadratcount(S, 20, 20))
  vL <- as.vector(quadratcount(L, 20, 20))
  
  gc_t <- geary.test(vT, listw = lw)$estimate[1]
  gc_s <- geary.test(vS, listw = lw)$estimate[1]
  gc_l <- geary.test(vL, listw = lw)$estimate[1]
  
  data.frame(GearyC.T = gc_t, GearyC.S = gc_s, GearyC.I = gc_l)
}

#-------------------------------------------------------------------------------
#' Bhattacharyya Coefficient Between Patterns
#'
#' Compute pairwise Bhattacharyya coefficient on flattened quadrat counts.
#'
#' @param T A `ppp` object of tumor points.
#' @param S A `ppp` object of stromal points.
#' @param L A `ppp` object of lymphocyte/immune points.
#' @return A 1-row `data.frame` with `BC.TS`, `BC.TI`, `BC.IS`.
#' @examples
#' bhattacharyya(tumor_ppp, stroma_ppp, lymph_ppp)
bhattacharyya <- function(T, S, L) {
  get_vec <- function(x) {
    v <- if (inherits(x, "ppp")) as.vector(quadratcount(x, 20, 20)) else as.numeric(x)
    as.numeric(v)
  }
  vT <- get_vec(T); vS <- get_vec(S); vL <- get_vec(L)
  
  bc_ts <- bhatt.coeff(vT, vS)
  bc_ti <- bhatt.coeff(vT, vL)
  bc_is <- bhatt.coeff(vL, vS)
  
  data.frame(BC.TS = bc_ts, BC.TI = bc_ti, BC.IS = bc_is)
}

#-------------------------------------------------------------------------------
#' Clark–Evans Nearest‐Neighbor Index
#'
#' Compute Clark–Evans R (Donnelly correction) for each pattern.
#'
#' @param T A `ppp` object of tumor points.
#' @param S A `ppp` object of stromal points.
#' @param L A `ppp` object of lymphocyte/immune points.
#' @return A 1-row `data.frame` with `CE.T`, `CE.S`, `CE.I`.
#' @examples
#' clark_evans(tumor_ppp, stroma_ppp, lymph_ppp)
clark_evans <- function(T, S, L) {
  
  CE.T <- clarkevans(T, correction = "Donnelly")
  CE.S <- clarkevans(S, correction = "Donnelly")
  CE.I <- clarkevans(L, correction = "Donnelly")
  
  data.frame(CE.T = CE.T, CE.S = CE.S, CE.I = CE.I)
}

#-------------------------------------------------------------------------------
#' Quadrat Count: VMR & Chi‐Square
#'
#' Compute Variance‐to‐Mean Ratio and Chi‐Square statistic on quadrat counts.
#'
#' @param T A `ppp` object of tumor points.
#' @param S A `ppp` object of stromal points.
#' @param L A `ppp` object of lymphocyte/immune points.
#' @return A 1-row `data.frame` with six columns:
#'   `T.VMR`, `T.ChiSq`, `S.VMR`, `S.ChiSq`, `L.VMR`, `L.ChiSq`.
#' @examples
#' quadrat_scalar_stats(tumor_ppp, stroma_ppp, lymph_ppp)
quadrat_scalar_stats <- function(T, S, L) {
  compute_stats <- function(ppp_obj) {
    counts <- as.vector(quadratcount(ppp_obj, 20, 20) %>% as.table())
    m <- mean(counts); v <- var(counts)
    vmr    <- v / m
    expect <- rep(m, length(counts))
    chi_sq <- sum((counts - expect)^2 / expect)
    c(VMR = vmr, ChiSq = chi_sq)
  }
  
  sT <- compute_stats(T); sS <- compute_stats(S); sL <- compute_stats(L)
  data.frame(
    T.VMR   = sT["VMR"],   T.ChiSq = sT["ChiSq"],
    S.VMR   = sS["VMR"],   S.ChiSq = sS["ChiSq"],
    L.VMR   = sL["VMR"],   L.ChiSq = sL["ChiSq"]
  )
}

#-------------------------------------------------------------------------------
#' Grid‐Wise Join Count Statistic
#'
#' Compute join counts for each pair of majority‐label cells in a 20×20 grid.
#'
#' @param T A `ppp` object of tumor points.
#' @param S A `ppp` object of stromal points.
#' @param L A `ppp` object of lymphocyte/immune points.
#' @return A 1-row `data.frame` with `JC.TS`, `JC.TI`, `JC.IS`.
#' @examples
#' join_count_stats(tumor_ppp, stroma_ppp, lymph_ppp)
join_count_stats <- function(T, S, L) {
  # label and merge
  marks(T) <- factor("T", levels = c("T","S","L"))
  marks(S) <- factor("S", levels = c("T","S","L"))
  marks(L) <- factor("L", levels = c("T","S","L"))
  all_ppp   <- superimpose(T, S, L, W = T$window)
  
  # grid‐wise majority label
  Q      <- quadratcount(all_ppp, 20, 20)
  labels <- apply(as.matrix(Q), 1:2, function(cell) {
    if (sum(cell)==0) return(NA)
    c("T","S","L")[which.max(cell)]
  }) 
  fac <- factor(as.vector(labels), levels = c("T","S","L"))
  
  # adjacency
  nb     <- cell2nb(20,20, type="queen")
  wmat   <- nb2mat(nb, style="B", zero.policy=TRUE)
  
  # count joins
  join_count <- function(lbl1, lbl2) {
    sum(sapply(seq_along(fac), function(i) {
      nei <- which(wmat[i,]==1)
      sum(fac[i]==lbl1 & fac[nei]==lbl2, na.rm=TRUE)
    }))
  }
  
  data.frame(
    JC.TS = join_count("T","S"),
    JC.TI = join_count("T","L"),
    JC.IS = join_count("L","S")
  )
}

#-------------------------------------------------------------------------------
#' Jaccard Index (mask out zero counts)
#'
#' Compute the generalized Jaccard similarity between each pair of point patterns
#' based on their 20×20 quadrat‐count vectors, ignoring cells where both counts are zero.
#'
#' J(A,B) = ⟨A,B⟩ / (‖A‖² + ‖B‖² − ⟨A,B⟩)
#'
#' @param T A `ppp` object of tumor points.
#' @param S A `ppp` object of stromal points.
#' @param L A `ppp` object of lymphocyte/immune points.
#' @return 1‐row `data.frame` with `Jaccard.TS`, `Jaccard.TI`, `Jaccard.IS`
#-------------------------------------------------------------------------------
jaccard_index <- function(T, S, L) {
  # helper: flatten 20×20 quadrat counts
  get_vec <- function(pp) {
    as.vector(as.table(quadratcount(pp, nx = 20, ny = 20)))
  }
  
  vT <- get_vec(T); vS <- get_vec(S); vL <- get_vec(L)
  
  calc_jaccard <- function(a, b) {
    # mask out cells where both are zero
    idx <- which(a > 0 | b > 0)
    a <- a[idx]; b <- b[idx]
    dot   <- sum(a * b)
    denom <- sum(a^2) + sum(b^2) - dot
    if (denom == 0) return(NA_real_)
    dot / denom
  }
  
  data.frame(
    Jaccard.TS = calc_jaccard(vT, vS),
    Jaccard.TI = calc_jaccard(vT, vL),
    Jaccard.IS = calc_jaccard(vL, vS)
  )
}

#-------------------------------------------------------------------------------
#' Dice–Sørensen Index (mask out zero counts)
#'
#' Compute the Dice (Sørensen) similarity between each pair of point patterns
#' based on their 20×20 quadrat‐count vectors, ignoring cells where both counts are zero.
#'
#' D(A,B) = 2⟨A,B⟩ / (‖A‖² + ‖B‖²)
#'
#' @param T A `ppp` object of tumor points.
#' @param S A `ppp` object of stromal points.
#' @param L A `ppp` object of lymphocyte/immune points.
#' @return 1‐row `data.frame` with `Dice.TS`, `Dice.TI`, `Dice.IS`
#-------------------------------------------------------------------------------
dice_sorensen_index <- function(T, S, L) {
  get_vec <- function(pp) {
    as.vector(as.table(quadratcount(pp, nx = 20, ny = 20)))
  }
  
  vT <- get_vec(T); vS <- get_vec(S); vL <- get_vec(L)
  
  calc_dice <- function(a, b) {
    idx <- which(a > 0 | b > 0)
    a <- a[idx]; b <- b[idx]
    num   <- 2 * sum(a * b)
    denom <- sum(a^2) + sum(b^2)
    if (denom == 0) return(NA_real_)
    num / denom
  }
  
  data.frame(
    Dice.TS = calc_dice(vT, vS),
    Dice.TI = calc_dice(vT, vL),
    Dice.IS = calc_dice(vL, vS)
  )
}

#-------------------------------------------------------------------------------
#' Cosine Similarity (mask out zero counts)
#'
#' Compute the cosine similarity between each pair of point patterns
#' based on their 20×20 quadrat‐count vectors, ignoring cells where both counts are zero.
#'
#' Cosine(A,B) = ⟨A,B⟩ / (‖A‖ · ‖B‖)
#'
#' @param T A `ppp` object of tumor points.
#' @param S A `ppp` object of stromal points.
#' @param L A `ppp` object of lymphocyte/immune points.
#' @return 1‐row `data.frame` with `Cosine.TS`, `Cosine.TI`, `Cosine.IS`
#-------------------------------------------------------------------------------
cosine_similarity <- function(T, S, L) {
  get_vec <- function(pp) {
    as.vector(as.table(quadratcount(pp, nx = 20, ny = 20)))
  }
  
  vT <- get_vec(T); vS <- get_vec(S); vL <- get_vec(L)
  
  calc_cosine <- function(a, b) {
    idx <- which(a > 0 | b > 0)
    a <- a[idx]; b <- b[idx]
    dot <- sum(a * b)
    na  <- sqrt(sum(a^2))
    nb  <- sqrt(sum(b^2))
    if (na == 0 || nb == 0) return(NA_real_)
    dot / (na * nb)
  }
  
  data.frame(
    Cosine.TS = calc_cosine(vT, vS),
    Cosine.TI = calc_cosine(vT, vL),
    Cosine.IS = calc_cosine(vL, vS)
  )
}











