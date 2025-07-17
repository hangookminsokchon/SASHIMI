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
plot(x = df$x, y = df$y, col = df$type, cex = 0.3, pch = 16)

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
K_function <- Ksingle(df$Tumor, df$Stromal, df$Immune)
plot(K_function$single.K.T)
plot(K_function$single.K.S)
plot(K_function$single.K.I)

# Single-type J-function: compute & plot the functional features
J_function <- Jsingle(df$Tumor, df$Stromal, df$Immune)
plot(J_function$single.J.T)
plot(J_function$single.J.S)
plot(J_function$single.J.I)

# Cross-type functions, pair correlation function, etc...
##  Other functional features work the same way
