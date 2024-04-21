# prepare section ----
# Load data
my_data <- load(file = "")

# Install and load required packages
# install.packages("dplyr")
# install.packages("tidyr")
# install.packages("gridExtra")
# install.packages("ggplot2")
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

# section 1 ----
# Question1.(i)(ii)

# Filter the required data
selected_drug_data <- drug_data %>%
  select(sample, auc)

selected_clinical_data <- clinical_data %>%
  select(sample, RUNX1)

# Based on the "sample column", merge two data frames
merged_clinical_drug_data <- merge(selected_drug_data, selected_clinical_data, 
                                   by = "sample")
# Divide into two groups based on the value of the RUNX1 

merged_clinical_drug_data_0 <- merged_clinical_drug_data %>%
  filter(RUNX1 == 0)

merged_clinical_drug_data_1 <- merged_clinical_drug_data %>%
  filter(RUNX1 == 1)

# Calculate the distribution and 95% confidence interval of AUC
summary_stats_0 <- summary(merged_clinical_drug_data_0$auc)
summary_stats_1 <- summary(merged_clinical_drug_data_1$auc)

ci_0 <- t.test(merged_clinical_drug_data_0$auc)$conf.int
ci_1 <- t.test(merged_clinical_drug_data_1$auc)$conf.int

# Output results
cat("Group with RUNX1=0 - AUC Distribution Summary:\n")
print(summary_stats_0)

cat("\n95% Confidence Interval for Group with RUNX1=0:\n")
print(ci_0)

cat("\nGroup with RUNX1=1 - AUC Distribution Summary:\n")
print(summary_stats_1)

cat("\n95% Confidence Interval for Group with RUNX1=1:\n")
print(ci_1)

# plot histogram
histogram_0 <- ggplot(merged_clinical_drug_data_0, aes(x = auc)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Distribution of AUC for RUNX1=0",
       x = "AUC",
       y = "Frequency")

histogram_1 <- ggplot(merged_clinical_drug_data_1, aes(x = auc)) +
  geom_histogram(binwidth = 1, fill = "green", color = "black", alpha = 0.7) +
  labs(title = "Distribution of AUC for RUNX1=1",
       x = "AUC",
       y = "Frequency")

# plot confidence intervals
ci_plot_0 <- ggplot(data = data.frame(ci = ci_0), aes(x = 1, ymin = ci[1], ymax = ci[2])) +
  geom_errorbar(y = mean(merged_clinical_drug_data_0$auc), color = "blue") +
  geom_point(y = mean(merged_clinical_drug_data_0$auc), color = "blue") +
  labs(title = "95% Confidence Interval for RUNX1=0",
       x = "",
       y = "AUC")

ci_plot_1 <- ggplot(data = data.frame(ci = ci_1), aes(x = 1, ymin = ci[1], ymax = ci[2])) +
  geom_errorbar(y = mean(merged_clinical_drug_data_1$auc), color = "green") +
  geom_point(y = mean(merged_clinical_drug_data_1$auc), color = "green") +
  labs(title = "95% Confidence Interval for RUNX1=1",
       x = "",
       y = "AUC")

# plot boxplot
boxplot_0 <- ggplot(merged_clinical_drug_data_0, aes(x = factor(1), y = auc, fill = factor(1))) +
  geom_boxplot(fill = "blue", alpha = 0.7) +
  labs(title = "Boxplot of AUC for RUNX1=0",
       x = "",
       y = "AUC")

boxplot_1 <- ggplot(merged_clinical_drug_data_1, aes(x = factor(1), y = auc, fill = factor(1))) +
  geom_boxplot(fill = "green", alpha = 0.7) +
  labs(title = "Boxplot of AUC for RUNX1=1",
       x = "",
       y = "AUC")

# Arrange graphs via grid.arrange function
grid.arrange(histogram_0, histogram_1, boxplot_0, boxplot_1, ci_plot_0, ci_plot_1, ncol = 2)


# Question1.(iii)
library(boot)
library(dplyr)

# Create a Bootstrap function, which is used to resample with replacement from the sample
bootstrap_function <- function(data, indices) {
  sampled_data <- data[indices, ]
  return(mean(sampled_data$auc))
}

# Set Bootstrap sampling times
num_bootstrap_samples <- 1000

# Bootstrap resampling
set.seed(123)  
bootstrap_with_RUNX1 <- boot(merged_clinical_drug_data_1, statistic = bootstrap_function, R = num_bootstrap_samples)

# Draw a histogram of Bootstrap results
histogram_Bootstrap_0 <- ggplot(data.frame(auc = bootstrap_with_RUNX1$t), aes(x = auc)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Bootstrap Distribution of AUC for RUNX1=1",
       x = "AUC",
       y = "Frequency")

grid.arrange(histogram_Bootstrap_0)

# section 2 ----
# Question2.(i)(ii)(iii)

# "selected_drug_data" contains two columns of data, "AUC" and "sample".
# Load necessary libraries

# install.packages("Hmisc")
#install.packages("htmltools", dependencies=TRUE)
library(tidyverse)
library("Hmisc")

# Transpose selected_drug_data
selected_drug_data_transpose<- t(selected_drug_data)

# Use sample as column name
colnames(selected_drug_data_transpose) <- selected_drug_data_transpose[1,] 

# Merge the two data sets to get data with "auc" and "gex"
combined_data_auc_gex <- rbind(selected_drug_data_transpose,gex)
combined_data_auc_gex <-  combined_data_auc_gex[-1, ]

# Transpose the “combined_data_auc_gex” data to match the requirements
combined_data_auc_gex__transpose <- t(combined_data_auc_gex)

# Calculate the correlation coefficient between auc and genes
regression_gene_auc <- rcorr(as.matrix(combined_data_auc_gex__transpose), type = c("pearson"))

# Extract the correlation coefficient and P value
cor_col <- regression_gene_auc$r[,1]
p_value_col <- regression_gene_auc$P[,1]
#head(p_value_col)
#head(cor_col)

# Store correlation coefficient and P value into dataframe
regression_result_df <- data.frame(
  gene = rownames(regression_gene_auc$r)[-1],  # The first row is AUC, so exclude
  cor = cor_col[-1],  # Exclude correlation coefficient with itself
  p_value = p_value_col[-1]  # Exclude and own P-value
)

# p_values is the original P value vector, which needs to be corrected using the BH method
adjusted_p_values <- p.adjust(regression_result_df$p_value, method = "BH")

# Add adjusted P-values to the results data frame
regression_result_df$adjusted_p_value <- adjusted_p_values

# Convert the adjusted P value to numeric type
regression_result_df$adjusted_p_value <- as.numeric(regression_result_df$adjusted_p_value)

# Print results
# head(regression_result_df)

# Histogram of P values before BH adjustment
hist(regression_result_df$p_value, main = "P Value Distribution (Unadjusted)",
     xlab = "P Value", ylab = "Frequency", col = "lightblue", border = "black")

# Histogram of P values after BH adjustment
hist(regression_result_df$adjusted_p_value, main = "Adjusted P Value Distribution (BH Method)",
     xlab = "Adjusted P Value", ylab = "Frequency", col = "lightgreen", border = "black")


# Screen for genes with high correlation
correlated04_result <- regression_result_df[regression_result_df$cor > 0.4, ]
head(correlated04_result)

# Filter rows with adjusted_p_value less than 0.05
significant005_results <- regression_result_df[regression_result_df$adjusted_p_value < 0.05, ]
head(significant005_results)


# Question2.(iv)

# Filter out the gene expression data of ST6GALNAC3
row_index <- rownames(gex) == "ST6GALNAC3"
gex_ST6GALNAC3 <- gex[row_index,]

# Create a new data frame with two columns: Sample and Value
gex_ST6GALNAC3 <- data.frame(
  sample = colnames(gex),
  ST6GALNAC3 = as.numeric(as.matrix(gex_ST6GALNAC3))
)

# Merge the gene expression data of ST6GALNAC3 with merged_clinical_drug_data
merged_ST6GALNAC3_auc_RUNX1 <- merge(merged_clinical_drug_data, gex_ST6GALNAC3, 
                                   by = "sample")

# Use the lm function to fit the linear regression model
model1 <- lm(auc ~ ST6GALNAC3, data = merged_ST6GALNAC3_auc_RUNX1)
model2 <- lm(auc ~ ST6GALNAC3 + RUNX1, data = merged_ST6GALNAC3_auc_RUNX1)
model3 <- lm(auc ~ RUNX1, data = merged_ST6GALNAC3_auc_RUNX1)

# View summary of regression model
summary(model1)
summary(model2)
summary(model3)



# section 3 ----
# Question3.(i)

# Calculate the distance matrix between samples
dist_matrix <- dist(t(gex), method = "euclidean")

# Perform hierarchical clustering
hclust_result <- hclust(dist_matrix, method = "ward.D2")

# Divide clusters based on hierarchical clustering results and select clusters at a specific level
clusters <- cutree(hclust_result, k = 5)

# Print the number of samples in each cluster
cluster_counts <- table(clusters)
for (i in seq_along(cluster_counts)) {
  cat("Cluster", i, ":", cluster_counts[i], "\n")
}

# Set graphics parameters and increase the size of the graph
par(mfrow = c(1, 1), mar = c(2,3,2,2) + 0.01)

# Draw hierarchical clustering dendrogram
p  <-  plot(hclust_result, main = "Hierarchical Clustering Dendrogram", cex = 0.2)

# Draw the boundaries of cluster divisions
rect.hclust(hclust_result, k = 5, border = 2:6)


# Question3.(ii)
#BiocManager::install("heatmaps")
library(heatmaps)
library("pheatmap")

# Transpose selected_clinical_data
selected_clinical_data_transpose<- t(selected_clinical_data)

# Use sample as the column name
colnames(selected_clinical_data_transpose) <- selected_clinical_data_transpose[1,] 

# Merge two data sets
combined_data_RUNX_gex <- rbind(selected_clinical_data_transpose,gex)
combined_data_RUNX_gex <-  combined_data_RUNX_gex[-1, ]

classification_info <- as.factor(combined_data_RUNX_gex[1, ])

# Creating a new dataframe for annotation
annotation_c <- data.frame(Classification = classification_info)

# Setting row names to match the column names of the original dataframe
rownames(annotation_c) <- colnames(combined_data_RUNX_gex)

pheatmap(gex, 
         cluster_rows = T,
         cluster_cols = T,
         annotation_col =annotation_c, # Sample classification data
         annotation_legend=TRUE, # Display sample classification
         show_rownames = F,
         show_colnames = F,
         scale = "none", # No normalization
         color =colorRampPalette(c("#8854d0", "#ffffff","#fa8231"))(100)
)


# section 4 ----
# Question4.(i)
# Directly use combined_data_auc_gex data in Question2

# Remove the first row (AUC value)
gene_expression_data <- combined_data_auc_gex[-1, ]

# Calculate the variance of genes
gene_variances <- apply(gene_expression_data, 1, var)

# Find the 100 genes with the highest variance
top_genes <- names(sort(gene_variances, decreasing = TRUE)[1:100])

# Select these 100 genes in the data frame
selected_genes_data <- combined_data_auc_gex[top_genes,]
selected_genes_data  <- rbind(selected_drug_data_transpose,selected_genes_data)
selected_genes_data <- selected_genes_data[-1,]

# Make sure selected_genes_data is a data frame
selected_genes_data <- t(selected_genes_data)
selected_genes_data <- as.data.frame(selected_genes_data)

# Use lapply to convert each column into a numeric type
#(Note: Not numeric type will cause a lot of trouble)
selected_genes_data[top_genes] <- lapply(selected_genes_data[top_genes], as.numeric)

# Define multiple linear regression model
lm_model <- lm(auc ~ ., data = selected_genes_data)

# Output model summary
summary(lm_model)

combined_data_auc_gex  <- as.data.frame(combined_data_auc_gex)
# Calculate the observed values and predicted values of sample BA2409
sample_index <- which(colnames(combined_data_auc_gex) == "BA2409")
observed_auc <- combined_data_auc_gex[sample_index,1]
predicted_auc <- predict(lm_model, newdata = selected_genes_data[sample_index, , drop = FALSE])

print(paste("Observed AUC:", observed_auc))
print(paste("predictedAUC:", predicted_auc))

# Get the observed values and predicted values of all samples
observed_auc_all <- combined_data_auc_gex[1, ]
predicted_auc_all <- predict(lm_model, newdata = selected_genes_data[,-1])  # 排除第一列（AUC）

# Convert observed values and predicted values to numeric types
observed_auc_all <- as.numeric(observed_auc_all)
predicted_auc_all <- as.numeric(predicted_auc_all)

# Draw the relationship between observed values and predicted values of all samples
plot(observed_auc_all, predicted_auc_all, main = "Observed value and predicted value relationship graph", 
     xlab = "Observed value", ylab = "predicted value", col = "blue", pch = 16)

# Add a diagonal line to represent perfect fit
abline(0, 1, col = "red")  

# Get the residuals of the model
residuals <- residuals(lm_model)

# Draw residual plot
plot(predicted_auc_all, residuals, main = "Residual Plot",
     xlab = "predicted AUC", ylab = "Residuals", col = "blue", pch = 16)

# Add a horizontal line to check the distribution of residuals
abline(h = 0, col = "red")  

# Draw a smooth curve and check the homogeneity of variances of the residuals
lines(lowess(predicted_auc_all, residuals), col = "green")

#Add legend
legend("topright", legend = c("Residuals", "Smoothed Line"), col = c("blue", "green"), pch = 16)


# Draw Q-Q plot
qqnorm(residuals)
qqline(residuals, col = "red")

# Add title and tags
xlabel <- "Theoretical Quantiles"
ylabel <- "Sample Quantiles"
xlab(xlabel)
ylab(ylabel)


# 4.(ii)

# Load required library
library(boot)

# Assuming loocv_model is your regression model
# predicted_auc_all is the predicted values
# observed_auc_all is the observed (actual) values
# selected_genes_data is your dataframe with AUC in the first column and other predictors in the remaining columns

# In-sample correlation coefficient
in_sample_cor <- cor(predicted_auc_all, observed_auc_all)

# Leave-One-Out Cross-Validation
# Initialize an empty vector to store cross-validated predictions
loocv_predictions <- numeric(length(observed_auc_all))


selected_genes_data <- as.data.frame(sapply(selected_genes_data, as.numeric))

str(selected_genes_data)

# Perform leave-one-out cross-validation
for (i in 1:length(observed_auc_all)) {
  # Exclude the i-th observation from the training set
  train_data <- selected_genes_data[-i, ]
  
  # Fit the model on the training data
  loocv_model <- lm(auc ~ ., data = train_data)
  
  # Predict on the i-th observation
  loocv_predictions[i] <- predict(loocv_model, newdata = selected_genes_data[i, , drop = FALSE])
}

# Calculate the correlation coefficient for leave-one-out cross-validation
loocv_cor <- cor(loocv_predictions, observed_auc_all)

# Display the results
cat("In-Sample Correlation Coefficient:", in_sample_cor, "\n")
cat("Leave-One-Out Cross-Validation Correlation Coefficient:", loocv_cor, "\n")


