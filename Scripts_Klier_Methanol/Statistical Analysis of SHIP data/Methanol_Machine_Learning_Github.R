#----------------------------------------------
## Install and load the required packages
#----------------------------------------------
library(glmnet)
library(imputeTS)
library(caret)
library("writexl")
library("readxl")
library(vip)  
library(boot)
library(mice)
library(stats)
library(ggplot2)
library(RColorBrewer)


#----------------------------------------------
## Import Data
#----------------------------------------------
setwd("C:/Users/mehrjerda/Desktop/Canada")
imputed_df <- read_excel("C:\\Users\\mehrjerda\\Desktop\\ML_Methanol\\Methanol_ML_Impu-Predictors_4p.xlsx")
summary(imputed_df)

#Rename the name of variables 
#imputed_df$resA_ln_meth_abso <- imputed_df

#names(imputed_df)[names(imputed_df) == "resA_ln_meth_abso"] <- "Methanol:2"
#imputed_df

#names(imputed_df)
#----------------------------------------------
# Handling missing values
#----------------------------------------------
# Impute missing values using mice
imputed_data <- mice(imputed_df, method = "pmm", m = 5, maxit = 50)

# Complete the imputed data
completed_data <- complete(imputed_data)

# Now you can use 'completed_data' for further analysis
data <-completed_data
#----------------------------------------------
##Predictors
#----------------------------------------------

Methanol_var <- completed_data[c(2:119)]
Methanol_len <- length(Methanol_var)
Methanol_var_names <- colnames(completed_data)[c(2:119)]

#----------------------------------------------
##Outcome_var
#----------------------------------------------
Outcome_var <-completed_data[c(1)]
Outcome_len <- length(Outcome_var)

Outcome_var_names <- colnames(completed_data)[c(1)]


#----------------------------------------------
## Separate the predictor variables (features) and the outcome variable
#----------------------------------------------

X <- data[, Methanol_var_names] %>%
  scale(center = TRUE, scale = FALSE) %>%
  as.matrix()

y <- data[[Outcome_var_names]]

#---------------------------------------------
##Export summary of the dataset
#---------------------------------------------

# Get the summary of the dataset using the 'summary' function
data_summary <- summary(data)

# Convert the summary results to a data frame
data_summary_df <- as.data.frame(data_summary)

# Specify the output file path and name
output_file <- "C:\\Users\\mehrjerda\\Desktop\\ML_Methanol\\summary_dataset.xlsx"

# Save the summary to an Excel file
write_xlsx(data_summary_df, output_file)

#-------------------------------------------------
##Distribution 
#--------------------------------------------
cor_matrix <- cor(data[, Methanol_var_names])
print(cor_matrix)

#---------------------------------
#Plot Heat Map to see correlation 
#----------------------------------
# Load the necessary library

# Create a heatmap for the correlation matrix
heatmap(cor_matrix, 
        col = colorRampPalette(c("#D7191C", "#FFFFBF", "#2C7BB6"))(100), 
        scale = "none", 
        symm = TRUE, 
        margins = c(20, 20))

##---------------------------------------
# Extract predictor variable names
Methanol_var_names <- colnames(completed_data)[c(2:119)]

# Extract outcome variable name
Outcome_var_name <- colnames(completed_data)[1]

# Separate the predictor variables (features) and the outcome variable
X <- completed_data[, Methanol_var_names] %>%
  scale(center = TRUE, scale = FALSE) %>%
  as.matrix()

y <- completed_data[[Outcome_var_name]]

# Calculate the correlation matrix
cor_matrix <- cor(X)

# Convert the correlation matrix to a data frame
cor_df <- as.data.frame(as.table(cor_matrix))

# Create a heatmap using ggplot2
ggplot(cor_df, aes(x = Var1, y = Var2, fill = Freq)) +
  geom_tile() +
  scale_fill_gradientn(colors = brewer.pal(n = 9, name = "RdYlBu")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right") +
  labs( x = "dietary intake frequency ", y = "dietary intake frequency ", fill = "Correlation")

#----------------------------------------------
## Split the data into train and test sets
#----------------------------------------------
set.seed(123)  # Set a seed for reproducibility
train_indices <- createDataPartition(as.vector(y), p = 0.8, list = FALSE)  # Convert 'y' to a vector
train_X <- X[train_indices, ]
train_y <- y[train_indices]
test_X <- X[-train_indices, ]
test_y <- y[-train_indices]


#-----------------------------------------------
# Cross Validation for obtaining best value of alpha and lambda for elasticnet regression 
#-----------------------------------------------
# Set up k-fold cross-validation with 10 folds
#control <- trainControl(
# method = "cv",
# number = 10,   # Number of folds (k)
# search = "random",
# verboseIter = TRUE
#)


# Set up Repeated CV 
control <- trainControl(
  method = "repeatedcv",
  number = 10,  # Increase the number of folds to your desired value
  repeats = 10,
  search = "random",
  verboseIter = TRUE
)


# Fit the elastic net model using the train data
elastic_model <- train(x = train_X,
                       y = train_y,
                       method = "glmnet",
                       preProcess = c("center", "scale"),
                       metric = "Rsquared",
                       tuneLength = 100,
                       trControl = control)

elastic_model


#--------------------------------------------------------------------
# Extract the cross-validation results from the elastic_model object
cv_results <- elastic_model$results

cv_results1 <-data.frame(cv_results)


write_xlsx(cv_results1, "C:\\Users\\mehrjerda\\Desktop\\ML_Methanol\\cv_results-EN.xlsx")


#--------------------------------------------------------------------
##Prediction vs. Observed  Plot Train and Test  
#--------------------------------------------------------------------

train_predicted <- predict(elastic_model, train_X)
test_predicted <- predict(elastic_model, test_X)
predicted_value <- predict(elastic_model, X)


#--------------------------------------------------------------------
## Extract Prediction Score 
#--------------------------------------------------------------------


# Step 2: Create a dataframe with observed values (y) and predicted values
result_Pred <- data.frame(Observed = y, Predicted = predicted_value)

# Step 2: Create a dataframe with observed values (y) and predicted values
write_xlsx(result_Pred, "C:\\Users\\mehrjerda\\Desktop\\ML_Methanol\\predicted_values-EN.xlsx")


predict(elastic_model, test_X)
#-------------------------------------------------------------------------------

# Create dataframes with observed and predicted values for train and test
train_obs_pred_df <- data.frame(Observed = train_y, Predicted = train_predicted)
test_obs_pred_df <- data.frame(Observed = test_y, Predicted = test_predicted)

# Plot observed vs. predicted for train data

p<- ggplot(train_obs_pred_df, aes(x =Predicted , y = Observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Methanol Concentration", y = "Predicted Methanol") +
  ggtitle("Train Data: Observed vs. Predicted") +
  theme_minimal()+ 
  theme(axis.title.x = element_text(size = 25), # Adjust the 'size' for x-axis label text
        axis.title.y = element_text(size = 25), # Adjust the 'size' for y-axis label text
        axis.text.x = element_text(size = 20), # Adjust the 'size' for x-axis tick labels
        axis.text.y = element_text(size = 20)) # Adjust the 'size' for y-axis tick labels

# Save the loading plots as image files
ggsave("C:\\Users\\mehrjerda\\Desktop\\ML_Methanol\\train.jpeg", p)

# Plot observed vs. predicted for test data
q <- ggplot(test_obs_pred_df, aes(x = Predicted, y = Observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Methanol Concentration", y = "Predicted Methanol") +
  ggtitle("Test Data: Observed vs. Predicted") +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 25), # Adjust the 'size' for x-axis label text
        axis.title.y = element_text(size = 25), # Adjust the 'size' for y-axis label text
        axis.text.x = element_text(size = 20), # Adjust the 'size' for x-axis tick labels
        axis.text.y = element_text(size = 20)) # Adjust the 'size' for y-axis tick labels
# Save the loading plots as image files
ggsave("C:\\Users\\mehrjerda\\Desktop\\ML_Methanol\\test.jpeg", q)



#-----------------------------------------------------
##Reposrt Measures based on outcome5
#-----------------------------------------------------
final_rmse <- elastic_model$results$RMSE[length(elastic_model$results$RMSE)]

final_rmse

#-----------------------------------------------
## Extract the RMSE, R-squared, and MAE from the elastic_model
#-----------------------------------------------
rmse <- elastic_model$results$RMSE
#best value of RMSE is smallest value
best_rmse <- min(rmse)


r_squared <- elastic_model$results$Rsquared
# Exclude NaN values and get the maximum R-squared value
max_r_squared <- max(r_squared[!is.nan(r_squared)])


mae <- elastic_model$results$MAE
#best value of MAE is smallest value
best_MAE<- min(mae)


#-----------------------------------------------
## Extract Alpha and Lambda from the elastic_model
#-----------------------------------------------
best_alpha <-elastic_model$bestTune$alpha 
best_lambda <-elastic_model$bestTune$lambda 

all_alpha <- cv_results$alpha
all_lambda <- cv_results$lambda
rmse <- cv_results$RMSE

mae<- cv_results$MAE

r_squared<- cv_results$Rsquared
#-----------------------------------------------
## Create a dataframe with the results
#-----------------------------------------------

result_df <- data.frame(alpha = all_alpha,
                        lambda = all_lambda,
                        rmse = rmse,
                        mae = mae,
                        r_squared = r_squared
)

colnames(result_df) <- c("Alpha", "Lambda", "RMSE", "MAE","R-squared" )
# Print the dataframe
print(result_df)

# Remove rows with NaN values from result_df
result_df_clean <- result_df[complete.cases(result_df), ]

# Extract the alpha, lambda, and R-squared values from the cleaned result_df data frame
alpha_values <- result_df_clean$Alpha
lambda_values <- result_df_clean$Lambda
r_squared_values <- result_df_clean$`R-squared`
# Find the maximum value of log(lambda)
max_log_lambda <- max(log(result_df$Lambda))

# Find the minimum value of log(lambda)
min_log_lambda <- min(log(result_df$Lambda))

# Print the results
print("Maximum value of log(lambda):")
print(max_log_lambda)

print("Minimum value of log(lambda):")
print(min_log_lambda)




# Create a data frame with the extracted values
plot_data <- data.frame(alpha = alpha_values, lambda = lambda_values, r_squared = r_squared_values)

# Create the scatter plot

# Create the scatter plot without alpha labels
ggplot(plot_data, aes(x = log(lambda), y = r_squared)) +
  geom_point(size = 3, color = ifelse(plot_data$alpha == 0, "yellow", "yellow")) +
  labs(x = "Log(Lambda)", y = "R-squared", title = "R-squared vs. Lambda by Alpha") +
  theme_minimal()



# Custom colors for Ridge and Lasso
ridge_color <- "#1f78b4"  # Blue color
lasso_color <- "#e31a1c"  # Red color

#-------------------------------------------------
# Create the scatter plot with custom colors
#-------------------------------------------------
#try to change the lable a-y axis 

# Define the theme with customized text sizes
custom_theme <- theme(
  axis.title.x = element_text(size = 25), # Adjust the 'size' for x-axis label text
  axis.title.y = element_text(size = 25), # Adjust the 'size' for y-axis label text
  axis.text.x = element_text(size = 20), # Adjust the 'size' for x-axis tick labels
  axis.text.y = element_text(size = 20) # Adjust the 'size' for y-axis tick labels
)

# Define your first plot with the custom theme
plambda <- ggplot(plot_data, aes(x = log(lambda), y = r_squared)) +
  geom_point(size = 3, color = lasso_color) +
  labs(x = "Log(Lambda)", y = "R-squared", title = "R-squared vs. Lambda") +theme_minimal()+
  custom_theme 

# Save the first plot
ggsave("C:\\Users\\mehrjerda\\Desktop\\ML_Methanol\\plambda.jpeg", plambda)

# Define your second plot with the custom theme
palpha <- ggplot(plot_data, aes(x = log(alpha), y = r_squared)) +
  geom_point(size = 3, color = ridge_color) +
  labs(x = "Log(Alpha)", y = "R-squared", title = "R-squared vs. Alpha") +
  theme_minimal()+ custom_theme

# Save the second plot
ggsave("C:\\Users\\mehrjerda\\Desktop\\ML_Methanol\\palpha.jpeg", palpha)


#-----------------------------------------------
##Export the dataframe as an Excel file
#-----------------------------------------------
write_xlsx(result_df, "C:\\Users\\mehrjerda\\Desktop\\ML_Methanol\\Reports_Measures_final-EN.xlsx")
#*---------------------------------------------------

# Load the required library
df<- read_excel("C:\\Users\\mehrjerda\\Desktop\\ML_Methanol\\Methanol-New-Annotation_last_4p.xlsx")
summary(df)

variable_names <- colnames(df)
# Extract the new names from the data frame and store them in the new_var_names vector


# Check if the length of new_var_names is 124
if (length(variable_names) != 118) {
  stop("The number of new variable names should be 124.")
}

# Now you have the new variable names in the new_var_names vector.
# You can use them to rename the variables in completed_data and X as shown in the previous response.
# Modify the column names of 'completed_data'
colnames(completed_data)[c(2:118)] <- variable_names

# Update the predictor variable names in 'X'
Methanol_var_names <- variable_names


# Fit the glmnet model using train data
fit <- glmnet(x = train_X, y = train_y, alpha = elastic_model$bestTune$alpha  , lambda = elastic_model$bestTune$lambda)
#--------------------------------------------
# Check for missing data in train_X (predictor variables)
#----------------------------------------------

missing_X <- colSums(is.na(train_X))

# Check for missing data in train_y (outcome variable)
missing_y <- sum(is.na(train_y))

# Display the results
print("Missing data in train_X:")
print(missing_X)

print("Missing data in train_y:")
print(missing_y)

#---------------------------------------------
# Extract coefficients and calculate confidence intervals
#----------------------------------------------

coefficients <- coef(fit)


coef_ci <- sapply(1:ncol(coefficients), function(i) {
  quantile(coefficients[, i], c(0.025, 0.975))
})

# Create a data frame with variable names, coefficients, and confidence intervals
coefficients_df <- data.frame(
  Feature = rownames(coefficients)[-1],  # Exclude the intercept (first coefficient)
  Coefficient = coefficients[-1],
  Lower_CI = coef_ci[1, ],
  Upper_CI = coef_ci[2, ]
)

# Print the data frame with coefficients and confidence intervals
print(coefficients_df)

#----------------------------------------------
# Assuming you have already loaded the necessary libraries (ggplot2)
# Create a plot with coefficients and confidence intervals
#----------------------------------------------

ggplot(coefficients_df, aes(x = Feature, y = Coefficient, ymin = Lower_CI, ymax = Upper_CI)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.5)) +
  labs(x = "Feature", y = "Coefficient", title = "Feature Coefficients with Confidence Intervals") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#----------------------------------------------
# Create a plot with coefficients and confidence intervals using geom_pointrange
#----------------------------------------------

ggplot(coefficients_df, aes(x = Feature, y = Coefficient, ymin = Lower_CI, ymax = Upper_CI)) +
  geom_pointrange() +
  labs(x = "Feature", y = "Coefficient", title = "Feature Coefficients with Confidence Intervals") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Alternatively, you can use geom_errorbar with wider error bars:
ggplot(coefficients_df, aes(x = Feature, y = Coefficient, ymin = Lower_CI, ymax = Upper_CI)) +
  geom_point() +
  geom_errorbar(width = 0.2) +
  labs(x = "Feature", y = "Coefficient", title = "Feature Coefficients with Confidence Intervals") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#---------------------------------------------------
# Step 1: Replace the names of features with the new names
coefficients_df$Feature <- Methanol_var_names[-1]
updated_feature_names <- c(Methanol_var_names[1], coefficients_df$Feature[-1])
# Step 2: Find the indices of old feature names that match the new feature names
match_indices <- match(updated_feature_names, coefficients_df$Feature)
coefficients_df$Feature <- updated_feature_names
# Step 3: Replace the old feature names with the new feature names
coefficients_df$Feature[!is.na(match_indices)] <- updated_feature_names[!is.na(match_indices)]

# Print the updated coefficients_df data frame
print(coefficients_df)


#---------------------------------------------------
# Specify the output file path and name
#---------------------------------------------------

output_file <- "C:\\Users\\mehrjerda\\Desktop\\ML_Methanol\\Coefficients_CI-EN.xlsx"

# Save the coefficients_df data frame to an Excel file
write_xlsx(coefficients_df, output_file)

#---------------------------------------------------
##Plot Importance of the features
#---------------------------------------------------
#---------------------------------------------------
#Feature-Importance based on filterd coefficent with zero values
#---------------------------------------------------
# Increase figure size
options(repr.plot.width = 10, repr.plot.height = 8)

# Plot the feature importance using ggplot2
feature_importance_plot <- ggplot(coefficients_df_filtered, aes(x = reorder(Feature, -Coefficient), y = Coefficient, fill = Coefficient > 0)) +
  geom_bar(stat = "identity", position = "identity", width = 0.9) +
  scale_fill_manual(values = my_colors, guide = FALSE) +
  labs(x = "Dietary Intake Frequency", y = "ElasticNet Coefficient", fill = "Coefficient > 0") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 1.2))  # Rotate x-axis labels

# Print the plot
print(feature_importance_plot)


# Save the second plot
ggsave("C:\\Users\\mehrjerda\\Desktop\\ML_Methanol\\feature_importance_plot.jpeg", feature_importance_plot)

#-----------------------------------------------------------
#Select some metabolite names to have better visualization 
#----------------------------------------------------------

# Adjust figure size
options(repr.plot.width = 10, repr.plot.height = 8)

# Subset the data to display only a subset of metabolite names
subset_df <- coefficients_df_filtered[1:20, ]  # Display only the top 20 metabolites

# Plot the feature importance using ggplot2 with horizontal bars
feature_importance_plot <- ggplot(subset_df, aes(y = reorder(Feature, Coefficient), x = Coefficient, fill = Coefficient > 0)) +
  geom_bar(stat = "identity", position = "identity", width = 0.9) +
  scale_fill_manual(values = my_colors, guide = FALSE) +
  labs(y = "Dietary Intake Frequency", x = "ElasticNet Coefficient", fill = "Coefficient > 0") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 1.2)) +  # Adjust text size for y-axis
  coord_flip()  # Flip the plot to have horizontal bars

# Print the plot
print(feature_importance_plot)

# Save the plot
ggsave("C:\\Users\\mehrjerda\\Desktop\\ML_Methanol\\feature_importance_plot_new.jpeg", feature_importance_plot, width = 10, height = 8)


