# Description: Compute a suite of spatial summary statistic functions for single
# - and cross-type analyses of Tumor, Stromal, and Immune point patterns.
#
# Dependencies: spatstat, spdep, dispRity, dplyr
#
# List of spatial summary statistics
# K-function family:
#   • Single-type K-function
#     - Neighborhood K-function
#     - Locally scaled K-function
#     - Directional K-function
#   • Cross-type K-function
#     - Local Cross-type K-function
# 
# G-function family:
#   • Single-type G-function
#   • Cross-type G-function
# 
# L-function family:
#   • Single-type L-function
#   • Cross-type L-function
# 
# J-function family:
#   • Single-type J-function
#   • Cross-type J-function
# 
# Pair correlation function family:
#   • Pair correlation function
#   • Multitype pair correlation function
# 
# Other functions:
#   • Marked Connection Function
#   • Multitype I-function
#-------------------------------------------------------------------------------

library(spatstat)    # core spatial functions

#-------------------------------------------------------------------------------
#' Single-type Ripley's K-function
#' 
#' Compute single-type K-function for each point pattern.
#' 
#' @param T,S,L `ppp` objects.
#' @param r Optional numeric vector of distances at which to evaluate.
#' @return A list with `Ksingle.T`, `Ksingle.S`, `Ksingle.I` (class `fv`).
K_single <- function(T, S, L, r = NULL) {
  list(
    single.K.T = Kest(T, r = r, correction = "Ripley"),
    single.K.S = Kest(S, r = r, correction = "Ripley"),
    single.K.I = Kest(L, r = r, correction = "Ripley")
  )
}

#-------------------------------------------------------------------------------
#' Neighborhood K-function (Local K)
#' 
#' Compute local (neighborhood) K-function for each point pattern.
#' 
#' @param T,S,L `ppp` objects.
#' @param r Numeric vector of distances.
#' @return A list with `localK.T`, `localK.S`, `localK.I` (class `fv` or matrix).
K_local <- function(T, S, L, r = NULL) {
  list(
    local.K.T = localK(T, r = r),
    local.K.S = localK(S, r = r),
    local.K.I = localK(L, r = r)
  )
}

#-------------------------------------------------------------------------------
#' Locally Scaled K-function
#' 
#' Compute scaled K(r)/(π r^2) for each point pattern.
#' 
#' @param T,S,L `ppp` objects.
#' @param r Optional distances.
#' @return A list with `Kscaled.T`, `Kscaled.S`, `Kscaled.I` (class `fv`).
K_scaled <- function(T, S, L, r = NULL) {
  scale_fun <- function(X) {
    K <- Kest(X, r = r, correction = "Ripley")
    K$scaled <- with(K, iso / (pi * r^2))
    K
  }
  list(
    scaled.K.T = scale_fun(T),
    scaled.K.S = scale_fun(S),
    scaled.K.I = scale_fun(L)
  )
}

#-------------------------------------------------------------------------------
#' Directional K-function (Sector K)
#' 
#' Compute sector (directional) K-function for each point pattern.
#' 
#' @param T,S,L `ppp` objects.
#' @param r Numeric vector of distances.
#' @param sector Numeric vector length 2: start and end angles in radians.
#' @return A list with `Ksector.T`, `Ksector.S`, `Ksector.I` (class `fv`).
K_sector <- function(T, S, L, r = NULL, sector = c(0, pi/2)) {
  list(
    sector.K.T = Ksector(T, r = r, sector = sector),
    sector.K.S = Ksector(S, r = r, sector = sector),
    sector.K.I = Ksector(L, r = r, sector = sector)
  )
}

#-------------------------------------------------------------------------------
#' Cross-type K-function
#' 
#' Compute cross-type K-function between each pair of patterns.
#' 
#' @param T,S,L `ppp` objects; marks will be set internally.
#' @param r Optional distances.
#' @return A list with `Kcross.TS`, `Kcross.TI`, `Kcross.IS` (class `fv`).
K_cross <- function(T, S, L, r = NULL) {
  marks(T) <- "T"; marks(S) <- "S"; marks(L) <- "L" # impose marks to point pattern
  all_ppp <- superimpose(T, S, L, W = T$window) # merges the point pattern into one
  all_ppp <- ppp(all_ppp$x, all_ppp$y, window = all_ppp$window,
                 marks = factor(all_ppp$marks)) 
  list(
    cross.K.TS = Kcross(all_ppp, i = "T", j = "S", r = r),
    cross.K.TI = Kcross(all_ppp, i = "T", j = "L", r = r),
    cross.K.IS = Kcross(all_ppp, i = "L", j = "S", r = r)
  )
}

#-------------------------------------------------------------------------------
#' Local Cross-type K-function
#' 
#' Compute local (neighborhood) cross-type K-function for each pair.
#' 
#' @param T,S,L `ppp` objects.
#' @param r Numeric vector of distances.
#' @return A list with `localKcross.TS`, `localKcross.TI`, `localKcross.IS`.
K_cross_local <- function(T, S, L, r = NULL) {
  marks(T) <- "T"; marks(S) <- "S"; marks(L) <- "L"
  all_ppp <- superimpose(T, S, L, W = T$window)
  all_ppp <- ppp(all_ppp$x, all_ppp$y, window = all_ppp$window,
                 marks = factor(all_ppp$marks))
  list(
    cross.local.K.TS = localKcross(all_ppp, i = "T", j = "S", r = r),
    cross.local.K.TI = localKcross(all_ppp, i = "T", j = "L", r = r),
    cross.local.K.IS = localKcross(all_ppp, i = "L", j = "S", r = r)
  )
}

#-------------------------------------------------------------------------------
#' Empty space F-Function
#' 
#' Compute F-function for each point pattern.
#' 
#' @param T,S,L `ppp` objects.
#' @param r Optional distances.
#' @return A list with `Fsingle.T`, `Fsingle.S`, `Fsingle.I` (class `fv`).
F_single <- function(T, S, L, r = NULL) {
  list(
    single.F.T = Fest(T, r = r),
    single.F.S = Fest(S, r = r),
    single.F.I = Fest(L, r = r)
  )
}

#-------------------------------------------------------------------------------
#' Nearest Neighbor G-Function
#' 
#' Compute G-function for each point pattern.
#' 
#' @param T,S,L `ppp` objects.
#' @param r Optional distances.
#' @return A list with `Gsingle.T`, `Gsingle.S`, `Gsingle.I` (class `fv`).
G_single <- function(T, S, L, r = NULL) {
  list(
    single.G.T = Gest(T, r = r),
    single.G.S = Gest(S, r = r),
    single.G.I = Gest(L, r = r)
  )
}

#-------------------------------------------------------------------------------
#' Cross-type G-function
#' 
#' Compute cross-type G-function between each pair.
#' 
#' @param T,S,L `ppp` objects.
#' @param r Optional distances.
#' @return A list with `Gcross.TS`, `Gcross.TI`, `Gcross.IS`.
G_cross <- function(T, S, L, r = NULL) {
  marks(T) <- "T"; marks(S) <- "S"; marks(L) <- "L"
  all_ppp <- superimpose(T, S, L, W = T$window)
  all_ppp <- ppp(all_ppp$x, all_ppp$y, window = all_ppp$window,
                 marks = factor(all_ppp$marks))
  list(
    cross.G.TS = Gcross(all_ppp, i = "T", j = "S", r = r),
    cross.G.TI = Gcross(all_ppp, i = "T", j = "L", r = r),
    cross.G.IS = Gcross(all_ppp, i = "L", j = "S", r = r)
  )
}

#-------------------------------------------------------------------------------
#' Single-type J-function
#' 
#' Compute Jest for each point pattern.
#' 
#' @param T,S,L `ppp` objects.
#' @param r Optional distances.
#' @return A list with `Jsingle.T`, `Jsingle.S`, `Jsingle.I` (class `fv`).
J_single <- function(T, S, L, r = NULL) {
  list(
    single.J.T = Jest(T, r = r),
    single.J.S = Jest(S, r = r),
    single.J.I = Jest(L, r = r)
  )
}

#-------------------------------------------------------------------------------
#' Cross-type J-function
#' 
#' Compute Jcross between each pair.
#' 
#' @param T,S,L `ppp` objects.
#' @param r Optional distances.
#' @return A list with `Jcross.TS`, `Jcross.TI`, `Jcross.IS`.
J_cross <- function(T, S, L, r = NULL) {
  marks(T) <- "T"; marks(S) <- "S"; marks(L) <- "L"
  all_ppp <- superimpose(T, S, L, W = T$window)
  all_ppp <- ppp(all_ppp$x, all_ppp$y, window = all_ppp$window,
                 marks = factor(all_ppp$marks))
  list(
    cross.J.TS = Jcross(all_ppp, i = "T", j = "S", r = r),
    cross.J.TI = Jcross(all_ppp, i = "T", j = "L", r = r),
    cross.J.IS = Jcross(all_ppp, i = "L", j = "S", r = r)
  )
}

#-------------------------------------------------------------------------------
#' Besag's L-function (single-type)
#' 
#' Compute Lest for each point pattern.
#' 
#' @param T,S,L `ppp` objects.
#' @param r Optional distances.
#' @return A list with `Lsingle.T`, `Lsingle.S`, `Lsingle.I` (class `fv`).
L_single <- function(T, S, L, r = NULL) {
  list(
    single.L.T = Lest(T, r = r),
    single.L.S = Lest(S, r = r),
    single.L.I = Lest(L, r = r)
  )
}

#-------------------------------------------------------------------------------
#' Cross-type L-function
#' 
#' Compute Lcross between each pair.
#' 
#' @param T,S,L `ppp` objects.
#' @param r Optional distances.
#' @return A list with `Lcross.TS`, `Lcross.TI`, `Lcross.IS`.
L_cross <- function(T, S, L, r = NULL) {
  marks(T) <- "T"; marks(S) <- "S"; marks(L) <- "L"
  all_ppp <- superimpose(T, S, L, W = T$window)
  all_ppp <- ppp(all_ppp$x, all_ppp$y, window = all_ppp$window,
                 marks = factor(all_ppp$marks))
  list(
    cross.L.TS = Lcross(all_ppp, i = "T", j = "S", r = r),
    cross.L.TI = Lcross(all_ppp, i = "T", j = "L", r = r),
    cross.L.IS = Lcross(all_ppp, i = "L", j = "S", r = r)
  )
}

#-------------------------------------------------------------------------------
#' Pair Correlation Function (single-type)
#' 
#' Compute pcf for each point pattern.
#' 
#' @param T,S,L `ppp` objects.
#' @param r Optional distances.
#' @return A list with `pcf.T`, `pcf.S`, `pcf.I` (class `fv`).
PairCorrelation <- function(T, S, L, r = NULL) {
  list(
    single.pcf.T = pcf(T, r = r),
    single.pcf.S = pcf(S, r = r),
    single.pcf.I = pcf(L, r = r)
  )
}

#-------------------------------------------------------------------------------
#' Cross-type Pair Correlation Function
#' 
#' Compute cross-type pcf between each pair.
#' 
#' @param T,S,L `ppp` objects.
#' @param r Optional distances.
#' @return A list with `pcfcross.TS`, `pcfcross.TI`, `pcfcross.IS`.
PairCorrelation_cross <- function(T, S, L, r = NULL) {
  marks(T) <- "T"; marks(S) <- "S"; marks(L) <- "L"
  all_ppp <- superimpose(T, S, L, W = T$window)
  all_ppp <- ppp(all_ppp$x, all_ppp$y, window = all_ppp$window,
                 marks = factor(all_ppp$marks))
  list(
    cross.pcf.TS = pcfcross(all_ppp, i = "T", j = "S", r = r),
    cross.pcf,TI = pcfcross(all_ppp, i = "T", j = "L", r = r),
    cross.pcf.IS = pcfcross(all_ppp, i = "L", j = "S", r = r)
  )
}

#-------------------------------------------------------------------------------
#' Cross-type I Function
#' 
#' Compute multitype I-function for each cross-type.
#' 
#' @param T,S,L `ppp` objects.
#' @param r Optional distances.
#' @return A list with `Icross.TS`, `Icross.TI`, `Icross.IS`.
I_cross <- function(T, S, L, r = NULL) {
  marks(T) <- "T"; marks(S) <- "S"; marks(L) <- "L"
  all_ppp <- superimpose(T, S, L, W = T$window)
  all_ppp <- ppp(all_ppp$x, all_ppp$y, window = all_ppp$window,
                 marks = factor(all_ppp$marks))
  list(
    cross.I.TS = Iest(all_ppp, i = "T", j = "S", r = r),
    cross.I.TI = Iest(all_ppp, i = "T", j = "L", r = r),
    cross.I.IS = Iest(all_ppp, i = "L", j = "S", r = r)
  )
}

#-------------------------------------------------------------------------------
#' Marked Connection Function (cross-type)
#' 
#' Compute markconnect for each cross-type.
#' 
#' @param T,S,L `ppp` objects.
#' @param r Numeric vector of distances.
#' @return A list with `markconnect.TS`, `markconnect.TI`, `markconnect.IS`.
MarkConnect_cross <- function(T, S, L, r = NULL) {
  marks(T) <- "T"; marks(S) <- "S"; marks(L) <- "L"
  all_ppp <- superimpose(T, S, L, W = T$window)
  all_ppp <- ppp(all_ppp$x, all_ppp$y, window = all_ppp$window,
                 marks = factor(all_ppp$marks))
  list(
    cross.markconnect.TS = markconnect(all_ppp, i = "T", j = "S", r = r),
    cross.markconnect.TI = markconnect(all_ppp, i = "T", j = "L", r = r),
    cross.markconnect.IS = markconnect(all_ppp, i = "L", j = "S", r = r)
  )
}

