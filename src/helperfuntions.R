# Description: Compute a suite of grid-based spatial statistics for Tumor, Stromal,
#   and Immune point patterns starting from raw coordinate+type data.
# Dependencies: spatstat, spdep, dispRity, dplyr
#-------------------------------------------------------------------------------

library(spatstat)        # ppp, owin, quadratcount, superimpose, clarkevans
library(spdep)           # cell2nb, nb2mat, nb2listw, Moran.I, lee.test, geary.test
library(dispRity)        # bhatt.coeff
library(dplyr)           # data manipulation

#-------------------------------------------------------------------------------
#' Normalize Coordinates to [0,1]
#'
#' @param df A data.frame with numeric columns `x`, `y`.
#' @return The same data.frame with `x`,`y` rescaled to [0,1].
#' @examples
#' normalize_coords(data.frame(x=runif(10,2,5), y=runif(10,-1,3)))
normalize_coords <- function(df) {
  df$x <- (df$x - min(df$x)) / (max(df$x) - min(df$x))
  df$y <- (df$y - min(df$y)) / (max(df$y) - min(df$y))
  df
}

#-------------------------------------------------------------------------------
#' Prepare Quadrat Counts for Each Cell Type
#'
#' @param df A data.frame with columns `x`,`y`,`type` (factor or character: "tumor","stroma","lymphocyte").
#' @return A named list of three numeric vectors (length 400) of 20×20 quadrat counts.
#' @examples
#' df <- data.frame(x=runif(50), y=runif(50), type=sample(c("tumor","stroma","lymphocyte"),50,T))
#' prepare_quadrat_counts(df)
prepare_quadrat_counts <- function(df, typelist = NA) {
  win <- owin(c(0,1), c(0,1))
  pp_all <- ppp(df$x, df$y, window=win)
  
  list(
    Tumor   = as.vector(quadratcount(pp_all[df$type=="tumor"],   nx=20, ny=20)),
    Stromal = as.vector(quadratcount(pp_all[df$type=="stroma"],  nx=20, ny=20)),
    Immune  = as.vector(quadratcount(pp_all[df$type=="lymphocyte"], nx=20, ny=20))
  )
}

#-------------------------------------------------------------------------------
#' Prepare ppp Objects for Each Cell Type
#'
#' @param df A data.frame with columns `x`,`y`,`type`.
#' @return A named list of three `ppp` objects for tumor, stroma, lymphocyte.
#' @examples
#' prepare_point_patterns(df)
prepare_point_patterns <- function(df, typelist = NA) {
  win <- owin(c(0,1), c(0,1))
  
  list(
    Tumor   = ppp(df$x[df$type=="tumor"],    df$y[df$type=="tumor"],    window=win),
    Stromal = ppp(df$x[df$type=="stroma"],   df$y[df$type=="stroma"],   window=win),
    Immune  = ppp(df$x[df$type=="lymphocyte"], df$y[df$type=="lymphocyte"], window=win)
  )
}

#-------------------------------------------------------------------------------
#' Master Pipeline: Calculate Areal Features
#'
#' @param raw_df A data.frame with numeric `x`,`y` and factor/character `type` ("tumor","stroma","lymphocyte").
#' @return A single-row data.frame combining Moran's I, Geary's C, Lee's L, Morisita-Horn,
#'   Bhattacharyya, Clark-Evans, quadrat VMR/Chi-Square, Jaccard, Dice, Cosine features.
#' @examples
#' calculate_areal_feature(df)
calculate_areal_feature <- function(raw_df) {
  # normalize if outside [0,1]
  raw_df <- normalize_coords(raw_df)
  raw_df$type <- factor(raw_df$type, levels=c("tumor","stroma","lymphocyte"))
  
  # prepare data structures
  ppps <- prepare_point_patterns(raw_df)
  
  # compute each feature set
  df_moran  <- moran_I(ppps$Tumor,  ppps$Stromal, ppps$Immune)
  df_geary  <- geary_C(ppps$Tumor,  ppps$Stromal, ppps$Immune)
  df_lee    <- lee_L(ppps$Tumor,    ppps$Stromal, ppps$Immune)
  df_mh     <- morisita_horn(ppps$Tumor, ppps$Stromal, ppps$Immune)
  df_bc     <- bhattacharyya(ppps$Tumor, ppps$Stromal, ppps$Immune)
  df_ce     <- clark_evans(ppps$Tumor, ppps$Stromal, ppps$Immune)
  df_qc     <- quadrat_scalar_stats(ppps$Tumor, ppps$Stromal, ppps$Immune)
  df_jacc   <- jaccard_index(ppps$Tumor, ppps$Stromal, ppps$Immune)
  df_dice   <- dice_sorensen_index(ppps$Tumor, ppps$Stromal, ppps$Immune)
  df_cos    <- cosine_similarity(ppps$Tumor, ppps$Stromal, ppps$Immune)
  
  # combine and return
  cbind(df_moran, df_geary, df_lee, df_mh, df_bc, df_ce,
        df_qc, df_jacc, df_dice, df_cos)
}
#-------------------------------------------------------------------------------


