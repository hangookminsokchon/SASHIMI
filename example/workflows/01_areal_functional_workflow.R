# This file serves as an example workflow of the framework: SASHIMI
# ------------------------------------------------------------------------------
# Import the raw data(point pattern data): same file as data/"example_point_patternA.csv
#                                                       data/"example_point_patternB.csv
df_A <- read.csv("example_point_patternA.csv")
df_B <- read.csv("example_point_patternB.csv")

# Set the column names
names(df_A) <- c("x", "y", "type")
names(df_B) <- c("x", "y", "type")
# ------------------------------------------------------------------------------
#' Ex)
#' Plot the raw df images
#' 
# Example image A
visual_point_pattern(df_A)

# Example image B
visual_point_pattern(df_B)
# ------------------------------------------------------------------------------
#' Ex)
#' Compute the areal features of 'df_A', 'df_B' images
#' areal features return a '1 x m dataframe' with scalar values
#' 
areal_features_A <- calculate_areal_feature(df_A)
areal_features_B <- calculate_areal_feature(df_B)

print(names(areal_features_A)) # list a suite of features
print(names(areal_features_B))

print(areal_features_A)        # print the result dataframe
print(areal_features_B)    
# ------------------------------------------------------------------------------
#' Ex)
#' Compute the functional features of 'df_A', 'df_B' images
#' functional features return a list of functional data, which can be ploted 
#' using plot().

# Preprocess the raw df to normalized point pattern data
df_A <- normalize_coords(df_A) # normalize the point pattern as max(x), max(y) = 1
df_B <- normalize_coords(df_B)

df_A <- prepare_point_patterns(df_A) # convert the csv data into a ppp object
df_B <- prepare_point_patterns(df_B)

# Single-type K-function: compute & plot the functional features
# plot_spatial_features() function generates a single visualization that includes all 6 functional graphs computed from df_A and df_B
# by changing the parameter 'feature_type=', the function computes different spatial summary statistics(functional) among the list below
plot_spatial_features(df_A, df_B, feature_type = "K_single")

# List of available functional features:
# K-function family:                      feature_type=:
#   • Single-type K-function               - "K_single"
#     - Neighborhood K-function            - "K_local"       
#     - Locally scaled K-function          - "K_scaled"
#     - Directional K-function             - "K_sector"
#   • Cross-type K-function                - "K_cross"
#     - Local Cross-type K-function        - "K_cross_local" 
# 
# G-function family:
#   • Single-type G-function               - "G_single"
#   • Cross-type G-function                - "G_cross"
# 
# L-function family:
#   • Single-type L-function               - "L_single"
#   • Cross-type L-function                - "L_cross"
# 
# J-function family:
#   • Single-type J-function               - "J_single"
#   • Cross-type J-function                - "J_cross"
# 
# Pair correlation function family:
#   • Pair correlation function            - "PairCorrelation"
#   • Multitype pair correlation function  - "PairCorrelation_cross"
# 
# Other functions:
#   • Marked Connection Function           - "MarkConnect_cross" 
#   • Multitype I-function                 - "I_cross"

## Other functional features can be plotted using plot_spatial_features() 
## Ex) plot_spatial_features(df_A, df_B, feature_type = "K_cross"),
##     plot_spatial_features(df_A, df_B, feature_type = "PairCorrelation"),
##     etc...

# ------------------------------------------------------------------------------
#' Ex)
#' Download the functional and areal features of 'df' image in .csv format
#' functional features return a list of functional data in 500 x n tabular data, whereas areal data returns a 1-row tabular data

# Download the K-fucntion feature data
# The same function introduced above, plot_spatial_features() also stores computed spatial summary statistics
# enabling direct feature download by accessing stored data.
extracted_features <- plot_spatial_features(df_A, df_B, feature_type = "K_single")

# Computed features from df_A
write.csv(extracted_features$dataset_A$single.K.T, file = 'appropriate output file path for areal features')
write.csv(extracted_features$dataset_A$single.K.S, file = 'appropriate output file path for areal features')
write.csv(extracted_features$dataset_A$single.K.I, file = 'appropriate output file path for areal features')

# Computed features from df_B
write.csv(extracted_features$dataset_B$single.K.T, file = 'appropriate output file path for areal features')
write.csv(extracted_features$dataset_B$single.K.S, file = 'appropriate output file path for areal features')
write.csv(extracted_features$dataset_B$single.K.I, file = 'appropriate output file path for areal features')


# Download the areal feature data
write.csv(areal_features_A, file = 'appropriate output file path for areal features') # Areal features are downloaded entirely 
write.csv(areal_features_B, file = 'appropriate output file path for areal features') # Areal features are downloaded entirely 




