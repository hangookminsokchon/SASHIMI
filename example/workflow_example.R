# Import the raw data(point pattern data)
df <- read.csv("downloads/breast_cancer_data.csv")
names(df) <- c("x", "y", "type")

# Run the full pipeline
features <- calculate_areal_feature(df)
print(features)
