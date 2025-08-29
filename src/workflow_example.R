# This file serves as an example workflow of the framework: SASHIMI
# ------------------------------------------------------------------------------
# Import the raw data(point pattern data): same file as data/ex_poinpattern.csv
df <- read.csv("/Users/yoolkyupark/Main/Qiwei_Lab/SASHIMI/data/breast_cancer_pointpattern.csv")

# Set the column names
names(df) <- c("x", "y", "type")

# ------------------------------------------------------------------------------
#' Ex)
#' Plot the raw df image
#' 
plot(x = df$x, y = df$y, col = factor(df$type), cex = 0.3, pch = 16)

# ------------------------------------------------------------------------------
#' Ex)
#' Compute the areal features of 'df' image
#' areal features return a '1 x m dataframe' with scalar values
#' 
areal_features <- calculate_areal_feature(df)
print(names(areal_features)) # list a suite of features
print(areal_features)        # print the result dataframe

# ------------------------------------------------------------------------------
#' Ex)
#' Compute the functional features of 'df' image
#' functional features return a list of functional data, which can be ploted 
#' using plot().

# Preprocess the raw df to normalized point pattern data
df <- normalize_coords(df) # normalize the point pattern as max(x), max(y) = 1
df <- prepare_point_patterns(df) # convert the csv data into a ppp object

# Single-type K-function: compute & plot the functional features
K_function <- K_single(df$Tumor, df$Stromal, df$Immune)
plot(K_function$single.K.T)
plot(K_function$single.K.S)
plot(K_function$single.K.I)

# Single-type J-function: compute & plot the functional features
J_function <- J_single(df$Tumor, df$Stromal, df$Immune)
plot(J_function$single.J.T)
plot(J_function$single.J.S)
plot(J_function$single.J.I)

# Cross-type functions, pair correlation function, etc...
##  Other functional features work the same way

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




