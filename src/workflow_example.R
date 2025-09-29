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
plot(x = df_A$x, y = df_A$y, col = factor(df_A$type), cex = 0.3, pch = 16)
plot(x = df_B$x, y = df_B$y, col = factor(df_B$type), cex = 0.3, pch = 16)
# ------------------------------------------------------------------------------
#' Ex)
#' Compute the areal features of 'df' image
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
#' Compute the functional features of 'df' image
#' functional features return a list of functional data, which can be ploted 
#' using plot().

# Preprocess the raw df to normalized point pattern data
df_A <- normalize_coords(df_A) # normalize the point pattern as max(x), max(y) = 1
df_B <- normalize_coords(df_B)

df_A <- prepare_point_patterns(df_A) # convert the csv data into a ppp object
df_B <- prepare_point_patterns(df_B)

# Single-type K-function: compute & plot the functional features
K_function_A <- K_single(df_A$Tumor, df_A$Stromal, df_A$Immune)
plot(K_function_A$single.K.T)
plot(K_function_A$single.K.S)
plot(K_function_A$single.K.I)

K_function_B <- K_single(df_B$Tumor, df_B$Stromal, df_B$Immune)
plot(K_function_B$single.K.T)
plot(K_function_B$single.K.S)
plot(K_function_B$single.K.I)

# Single-type J-function: compute & plot the functional features
J_function_A <- J_single(df_A$Tumor, df_A$Stromal, df_A$Immune)
plot(J_function_A$single.J.T)
plot(J_function_A$single.J.S)
plot(J_function_A$single.J.I)

J_function_B <- J_single(df_B$Tumor, df_B$Stromal, df_B$Immune)
plot(J_function_B$single.J.T)
plot(J_function_B$single.J.S)
plot(J_function_B$single.J.I)

# Cross-type functions, pair correlation function, etc...
## Other functional features work the same way

# ------------------------------------------------------------------------------
#' Ex)
#' Download the functional and areal features of 'df' image in .csv format
#' functional features return a list of functional data in 500 x 5 tabular data, whereas areal data returns a 1-row tabular data

# Download the K-fucntion feature data
write.csv(K_function$single.K.T, file = 'appropriate output file path for single.K.T')
write.csv(K_function$single.K.S, file = 'appropriate output file path for single.K.S')
write.csv(K_function$single.K.I, file = 'appropriate output file path for single.K.I')

# Download the areal feature data
write.csv(areal_features, file = 'appropriate output file path for areal features') # Areal features are downloaded entirely 




