# Installing the packages
packages <- c("tidyverse", "ggplot2", "pheatmap", "ggfortify", 
              "caret", "randomForest", "kernlab", "kknn", "xgboost",
              "nnet", "e1071") 
install.packages(packages)

# Loading the Libraries
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(ggfortify)
library(caret)
library(randomForest)
library(kernlab)
library(kknn)
library(xgboost)
library(nnet)
library(e1071)

#STEP 1 END

#STEP 2

data <- read.csv("metabolite_data.csv", header = TRUE)
view(data)


#separate numeric and non-numeric columns for further processing
numeric_data <- data %>% select(where(is.numeric))
non_numeric_data <- data %>% select(where(negate(is.numeric)))


#STEP 3: DATA CLEANING FFORM
#Identify and print the number of missing values in each numeric column
missing_summary <- sapply(numeric_data, function(x) sum(is.na(x)))
missing_summary

#Replace missing values in numeric data
numeric_data <- numeric_data %>%
  mutate(across(everything(), ~ifelse(is.na(.), median(.,na.rm =TRUE), .)))

  

#Identify and handle negative values in numeric data
if(any(numeric_data < 0, na. = TRUE)) {
  warning("Negative values detected setting to NA.")
  numeric_data[numeric_data<0] <- NA
  numeric_data <- numeric_data%>%
    mutate(across(everything(), ~ifelse(is.na(.), median(.,na.rm =TRUE), .)))
} 

#step 4: Summary statistics

#compute summary statistics(mean, median, sd) for numeric data

summary_stats <- numeric_data %>%
  summarise(across(everything(),list(mean=mean, median = median, sd=sd)))
summary_stats
view(summary_stats)

#step 5 : visualization

#Generate heatmap of non numeric data
pheatmap(as.matrix(numeric_data), scale="row",
         main = "Heatmap of Metabolites", cluster_rows = TRUE, cluster_cols = TRUE)

#Perform PCA on numeric data and visualize the results
pca <- prcomp(numeric_data, center = TRUE, scale = TRUE)
autoplot(pca, main = "PCA plot")

# Convert data from wide to long format
numeric_long <- pivot_longer(numeric_data, cols = everything(), names_to = "metabolite", values_to = "values")

# Create box plot
ggplot(numeric_long, aes(x = metabolite, y = values)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Boxplot of Metabolites")

#step 6: statistical analysis

# add a stimulated group column for demonstration
set.seed(123)
data$Group <- sample(c("A", "B"), nrow(data), replace = TRUE)

#Perform Anova for each metabolite to test for differences across groups
anova_results <- sapply(names(numeric_data), function(metabolite){
  formula <- as.formula(paste(metabolite, "~Group"))
  summary(aov(formula, data = cbind(numeric_data, Group = data$Group)))
})
anova_results

# write Anova results to a csv file for a reference
write.csv(anova_results, "anova_results.csv", row.names = TRUE)

# Step 7: Machine Learning: Comparing Algorithms
# --------------------------------------------
# Check if the Group column exists for classification tasks
if ("Group" %in% colnames(data)) {
  # Ensure the Group column is treated as a categorical factor
  data$Group <- as.factor(data$Group)
  # Prepare data for modeling by combining numeric features and the Group column
  modeling_data <- data %>% select(where(is.numeric), Group)
  # Split data into training (80%) and testing (20%) sets
  set.seed(123)
  trainIndex <- createDataPartition(modeling_data$Group, p = 0.8, list = FALSE)
  trainData <- modeling_data[trainIndex, ]
  testData <- modeling_data[-trainIndex, ]
  # Specify machine learning algorithms to compare
  algorithms <- c("rf", "svmLinear", "knn", "nnet", "xgbTree")
  # Train and evaluate models for each algorithm
  models <- list()
  for (algo in algorithms) {
    set.seed(123)
    if (algo == "xgbTree") {
      # Special handling for xgboost model
      dtrain <- xgb.DMatrix(data = as.matrix(trainData %>% select(-Group)), label = as.numeric(trainData$Group) - 1)
      param <- list(objective = "binary:logistic", eval_metric = "error")
      models[[algo]] <- xgb.train(params = param, data = dtrain, nrounds = 100)
    } else {
      # Train other models using caret
      model <- train(
        Group ~ ., 
        data = trainData, 
        method = algo, 
        trControl = trainControl(method = "cv", number = 10),
        tuneLength = 5
      )
      models[[algo]] <- model
    }
    print(paste("Trained", algo))
  }
  # Evaluate and identify the best-performing model
  best_model_name <- NA
  highest_accuracy <- -Inf
  for (algo in names(models)) {
    if (algo == "xgbTree") {
      # Evaluate xgboost model
      dtest <- xgb.DMatrix(data = as.matrix(testData %>% select(-Group)))
      xgb_predictions <- predict(models[[algo]], newdata = dtest)
      xgb_predictions <- factor(ifelse(xgb_predictions > 0.5, levels(testData$Group)[2], levels(testData$Group)[1]))
      xgb_confusion <- confusionMatrix(xgb_predictions, testData$Group)
      xgb_accuracy <- xgb_confusion$overall["Accuracy"]
      if (xgb_accuracy > highest_accuracy) {
        highest_accuracy <- xgb_accuracy
        best_model_name <- algo
      }
    } else {
      # Evaluate models trained using caret
      caret_accuracy <- max(models[[algo]]$results$Accuracy, na.rm = TRUE)
      if (caret_accuracy > highest_accuracy) {
        highest_accuracy <- caret_accuracy
        best_model_name <- algo
      }
    }
  }
  if (!is.na(best_model_name)) {
    best_model <- models[[best_model_name]]
    print(paste("Best model is:", best_model_name))
    # Generate predictions using the best model
    if (best_model_name == "xgbTree") {
      predictions <- xgb_predictions
    } else {
      predictions <- predict(best_model, newdata = testData)
    }
    # Evaluate the best model using a confusion matrix
    confusion <- confusionMatrix(predictions, testData$Group)
    print(confusion)
  } else {
    stop("No valid best model could be identified.")
  }
} else {
  message("No 'Group' column found in the data. Skipping machine learning analysis.")
}

# --------------------------------------------
# Step 8: Model Testing Section
# --------------------------------------------
# Generate synthetic data for testing
set.seed(456)
synthetic_data <- tibble(
  Feature1 = rnorm(200, mean = 5, sd = 1),
  Feature2 = rnorm(200, mean = 10, sd = 2),
  Group = factor(sample(c("A", "B"), 200, replace = TRUE))
)

# Display the first few rows of synthetic data
print(head(synthetic_data))

# Split synthetic data into training (80%) and testing (20%) sets
synthetic_trainIndex <- createDataPartition(synthetic_data$Group, p = 0.8, list = FALSE)
synthetic_trainData <- synthetic_data[synthetic_trainIndex, ]
synthetic_testData <- synthetic_data[-synthetic_trainIndex, ]

# Train a random forest model on synthetic data
synthetic_model <- train(
  Group ~ ., 
  data = synthetic_trainData, 
  method = "knn", 
  trControl = trainControl(method = "cv", number = 10),
  tuneLength = 5
)

# Test the model on synthetic testing data
synthetic_predictions <- predict(synthetic_model, newdata = synthetic_testData)

# Evaluate the model performance using a confusion matrix
synthetic_confusion <- confusionMatrix(synthetic_predictions, synthetic_testData$Group)
print(synthetic_confusion)

# Interpret model performance
cat("Model trained on synthetic data shows an accuracy of", 
    synthetic_confusion$overall["Accuracy"], 
    "and a kappa of", 
    synthetic_confusion$overall["Kappa"], 
    "\nThis indicates the model's ability to generalize the classification.")


