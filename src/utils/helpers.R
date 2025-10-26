# Description: Compute a suite of grid-based spatial statistics for Tumor, Stromal,
#   and Immune point patterns starting from raw coordinate+type data.
# Dependencies: spatstat, spdep, dispRity, dplyr
#-------------------------------------------------------------------------------

library(spatstat)        # ppp, owin, quadratcount, superimpose, clarkevans
library(spdep)           # cell2nb, nb2mat, nb2listw, Moran.I, lee.test, geary.test
library(dispRity)        # bhatt.coeff
library(dplyr)           # data manipulation
#-------------------------------------------------------------------------------
#' Standardize Cell Type Names
#'
#' @param df A data.frame with a column `type` containing cell type labels.
#' @return The same data.frame with standardized `type` values: 
#'         'tumor', 'stroma', or 'lymphocyte'.
#' @details Performs case-insensitive matching to handle variations like:
#'   - Stromal, Stroma, stroma, stromal → 'stroma'
#'   - Tumor, tumor cell, Tumor cell → 'tumor'  
#'   - lymphocyte, immune, immune cell → 'lymphocyte'
#' @examples
#' df <- data.frame(x=1:5, y=1:5, type=c("Tumor", "stroma", "Immune cell", "TUMOR", "Stromal"))
#' standardize_cell_types(df)
standardize_cell_types <- function(df) {
  # Create a working copy to avoid modifying original
  df$type <- as.character(df$type)
  
  # Case-insensitive pattern matching with grepl
  # Order matters: check most specific patterns first
  
  # Stromal variations
  stroma_pattern <- grepl("stroma", df$type, ignore.case = TRUE)
  df$type[stroma_pattern] <- "stroma"
  
  # Tumor variations  
  tumor_pattern <- grepl("tumor", df$type, ignore.case = TRUE)
  df$type[tumor_pattern] <- "tumor"
  
  # Immune/lymphocyte variations
  immune_pattern <- grepl("lymphocyte|immune", df$type, ignore.case = TRUE)
  df$type[immune_pattern] <- "lymphocyte"
  
  # Check for unmatched types
  valid_types <- c("tumor", "stroma", "lymphocyte")
  unmatched <- !df$type %in% valid_types
  
  if (any(unmatched)) {
    warning(sprintf("Found %d cells with unrecognized types: %s\nThese will be excluded from analysis.",
                    sum(unmatched),
                    paste(unique(df$type[unmatched]), collapse=", ")))
    df <- df[!unmatched, ]
  }
  
  return(df)
}

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
#' Prepare ppp Objects for Each Cell Type
#'
#' @param df A data.frame with columns `x`,`y`,`type`.
#' @return A named list of three `ppp` objects for tumor, stroma, lymphocyte.
#' @examples
#' prepare_point_patterns(df)
prepare_point_patterns <- function(df, typelist = NA) {
  win <- owin(c(0,1), c(0,1))

  df <- standardize_cell_types(df)
  
  list(
    Tumor   = ppp(df$x[df$type=="tumor"],    df$y[df$type=="tumor"],    window=win),
    Stromal = ppp(df$x[df$type=="stroma"],   df$y[df$type=="stroma"],   window=win),
    Immune  = ppp(df$x[df$type=="lymphocyte"], df$y[df$type=="lymphocyte"], window=win)
  )
}
#-------------------------------------------------------------------------------
#' Visualize Point Pattern Data
#'
#' Creates a spatial point pattern visualization with points colored by cell type
#' and a legend on the right side.
#'
#' @param df A data.frame with columns `x`, `y`, `type` where:
#'   - `x`: x-coordinates of points
#'   - `y`: y-coordinates of points
#'   - `type`: cell type labels (factor or character)
#' @param point_size Numeric, size of points (default = 0.3)
#' @param output_file Optional character string for output PNG filename. 
#'   If NULL, plots to current device.
#' @param width Numeric, width of output image in inches (default = 8)
#' @param height Numeric, height of output image in inches (default = 6)
#' @param dpi Numeric, resolution for PNG output (default = 300)
#' 
#' @return Invisibly returns the factor levels of cell types.
#' 
#' @examples
#' # Plot to screen
#' visual_point_pattern(df_A)
#' 
#' # Save to file
#' visual_point_pattern(df_A, output_file = "tissue_A.png")
#' 
#' # Custom point size
#' visual_point_pattern(df_B, point_size = 0.5)
#'
#' @export
visual_point_pattern <- function(df, 
                                  point_size = 0.3,
                                  output_file = NULL,
                                  width = 8,
                                  height = 6,
                                  dpi = 300) {
  
  # Input validation
  if (!is.data.frame(df)) {
    stop("Input must be a data.frame")
  }
  
  required_cols <- c("x", "y", "type")
  if (!all(required_cols %in% names(df))) {
    stop("Data frame must contain columns: x, y, type")
  }
  
  if (nrow(df) == 0) {
    stop("Data frame is empty")
  }
  
  # Convert type to factor
  cell_types <- factor(df$type)
  n_types <- nlevels(cell_types)
  
  # Open PNG device if output file specified
  if (!is.null(output_file)) {
    png(output_file, width = width, height = height, units = "in", res = dpi)
  }
  
  # Set up layout: left for points (75%), right for legend (25%)
  layout(matrix(c(1, 2), nrow = 1), widths = c(0.75, 0.25))
  
  # Plot points
  par(mar = c(0, 0, 0, 0), bg = "white")
  plot(x = df$x, y = df$y,
       col = cell_types,  # Automatically assigns the color for each cell type
       cex = point_size, 
       pch = 16, 
       asp = 1,
       ann = FALSE, 
       axes = FALSE)
  
  # Plot legend
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", 
         legend = levels(cell_types),
         col = 1:n_types,  # Basic color palette
         pch = 16, 
         bty = "n", 
         title = "Cell Type",
         cex = 1.0,
         y.intersp = 1.2)
  
  # Close device if file output
  if (!is.null(output_file)) {
    dev.off()
    message(paste("Plot saved to:", output_file))
  }
  
  # Reset layout
  layout(1)
  
  # Return cell type levels invisibly
  invisible(levels(cell_types))
}
#-------------------------------------------------------------------------------
#' Plot spatial features for two datasets
#' 
#' Compare spatial summary statistics between two point pattern datasets
#' by plotting selected features side by side.
#' 
#' @param df_A First dataset (list with Tumor, Stromal, Immune ppp objects)
#' @param df_B Second dataset (list with Tumor, Stromal, Immune ppp objects)
#' @param feature_type Character string specifying which feature to compute.
#'   Options: "K_single", "K_cross", "K_local", "K_scaled", "K_sector",
#'           "F_single", "G_single", "G_cross", "J_single", "J_cross",
#'           "L_single", "L_cross", "PairCorrelation", "PairCorrelation_cross",
#'           "I_cross", "MarkConnect_cross", "K_cross_local"
#' @param r Optional numeric vector of distances
#' @param ... Additional parameters for specific functions (e.g., sector for K_sector)
#' @return Invisible NULL (plots are generated as side effect)
plot_spatial_features <- function(df_A, df_B, feature_type, r = NULL) {
  # Source objects(functions) from src/features/functional.R
  # if doesn't work, use absolute path instead
  # ex) source("C:/Users/YourName/Documents/SASHIMI/src/features/functional.R")
  source("src/features/functional.R")
  
  # Get the appropriate function
  feature_func <- get(feature_type)
  
  # Compute features for both datasets
  features_A <- feature_func(df_A$Tumor, df_A$Stromal, df_A$Immune, r = NULL)
  features_B <- feature_func(df_B$Tumor, df_B$Stromal, df_B$Immune, r = NULL)
  
  # Set up plotting layout - 2 rows (A and B), 3 columns (for each cell type or pair)
  par(mfrow = c(2, 3))
  
  # Plot each feature
  for (i in 1:length(features_A)) {
    plot(features_A[[i]], main = paste("Dataset A:", names(features_A)[i]))
  }
  for (i in 1:length(features_B)) {
    plot(features_B[[i]], main = paste("Dataset B:", names(features_B)[i]))
  }
  
  # Reset plotting parameters
  par(mfrow = c(1, 1))

  # Store computed features
  return(list(
    dataset_A = features_A,
    dataset_B = features_B
  ))
  
  invisible(NULL)
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


