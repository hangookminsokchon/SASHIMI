####################################################
########## Point Pattern Data Preparation ##########
####################################################


###### Normalize Coordinates to [0, 1] if Needed ######
normalize_coords <- function(df) {
  df$x <- (df$x - min(df$x)) / (max(df$x) - min(df$x))
  df$y <- (df$y - min(df$y)) / (max(df$y) - min(df$y))
  return(df)
}


###### Prepare Vectorized Quadrat Counts ######
prepare_quadrat_counts <- function(df) {
  win <- owin(c(0, 1), c(0, 1))
  pp_all <- ppp(x = df$x, y = df$y, window = win)
  
  q_tumor   <- as.vector(quadratcount(pp_all[df$type == "tumor"],   nx = 20, ny = 20))
  q_stromal <- as.vector(quadratcount(pp_all[df$type == "stroma"], nx = 20, ny = 20))
  q_immune  <- as.vector(quadratcount(pp_all[df$type == "lymphocyte"],  nx = 20, ny = 20))
  
  return(list(
    Tumor   = q_tumor,
    Stromal = q_stromal,
    Immune  = q_immune
  ))
}


###### Prepare ppp Objects for Each Cell Type ######
prepare_point_patterns <- function(df) {
  win <- owin(c(0, 1), c(0, 1))
  
  pp_tumor   <- ppp(df$x[df$type == "tumor"],   df$y[df$type == "tumor"],   window = win)
  pp_stromal <- ppp(df$x[df$type == "stroma"], df$y[df$type == "stroma"], window = win)
  pp_immune  <- ppp(df$x[df$type == "lymphocyte"],  df$y[df$type == "lymphocyte"],  window = win)
  
  return(list(
    Tumor   = pp_tumor,
    Stromal = pp_stromal,
    Immune  = pp_immune
  ))
}


###### Master Pipeline: Calculate Areal Features ######
calculate_areal_feature <- function(raw_df) {
  # Normalize coordinates if not in [0, 1]
  if (max(raw_df$x) > 1 || max(raw_df$y) > 1 || min(raw_df$x) < 0 || min(raw_df$y) < 0) {
    raw_df <- normalize_coords(raw_df)
  }
  
  # Ensure type is a factor with correct levels
  raw_df$type <- factor(raw_df$type)
  
  # Extract data representations
  quads <- prepare_quadrat_counts(raw_df)
  ppps  <- prepare_point_patterns(raw_df)
  
  # Compute areal features
  df_moran <- moran_I(quads$Tumor, quads$Stromal, quads$Immune)
  df_geary <- geary_C(quads$Tumor, quads$Stromal, quads$Immune)
  df_lee   <- lee_L(quads$Tumor, quads$Stromal, quads$Immune)
  df_mh    <- morisita_horn(quads$Tumor, quads$Stromal, quads$Immune)
  df_bc    <- bhattacharyya(quads$Tumor, quads$Stromal, quads$Immune)
  df_ce    <- clark_evans(raw_df)
  
  # Return combined feature data frame
  return(cbind(df_moran, df_geary, df_lee, df_mh, df_bc, df_ce))
}

