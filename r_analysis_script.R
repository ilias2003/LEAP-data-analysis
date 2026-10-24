# ==============================================================================
# LEAP Program - Scenario 1: Understanding Diabetes Progression
# Complete R Script for Data Cleaning, EDA, Modeling, and Diagnostics
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Helper Functions and Setup
# ------------------------------------------------------------------------------

# Custom function for descriptive statistics
describe_numeric <- function(x) {
  x_clean <- x[!is.na(x)]
  c(
    Mean   = mean(x_clean),
    SD     = sd(x_clean),
    Median = median(x_clean),
    IQR    = IQR(x_clean),
    Min    = min(x_clean),
    Max    = max(x_clean)
  )
}

# Custom function to calculate Variance Inflation Factor (VIF)
# (Avoids dependency on external packages like `car`)
manual_vif <- function(model) {
  X <- model.matrix(model)
  if ("(Intercept)" %in% colnames(X)) {
    X <- X[, -1, drop = FALSE]
  }
  vif_vals <- numeric(ncol(X))
  names(vif_vals) <- colnames(X)
  
  for (i in seq_len(ncol(X))) {
    y_temp <- X[, i]
    X_temp <- X[, -i, drop = FALSE]
    if (ncol(X_temp) > 0) {
      r_sq <- summary(lm(y_temp ~ X_temp))$r.squared
      vif_vals[i] <- 1 / (1 - r_sq)
    } else {
      vif_vals[i] <- 1
    }
  }
  return(vif_vals)
}

# ------------------------------------------------------------------------------
# 1. Question 1: Understand and Prepare the Data
# ------------------------------------------------------------------------------

# Load the dataset (download if not present locally)
file_path <- "diabetes.tab.txt"
if (!file.exists(file_path)) {
  download.file("https://www2.stat.duke.edu/~rcs48/data/diabetes.tab.txt", destfile = file_path)
}

diabetes_raw <- read.table(file_path, header = TRUE, sep = "\t")

# Check basic structure and dimensions
cat("=== Data Structure ===\n")
str(diabetes_raw)
cat("\nDimensions (Rows x Columns):", dim(diabetes_raw), "\n")

# Check for missing values and duplicate records
cat("\nMissing values per variable:\n")
print(colSums(is.na(diabetes_raw)))

cat("\nNumber of duplicate rows:", sum(duplicated(diabetes_raw)), "\n")

# Data cleaning & re-encoding
diabetes <- diabetes_raw
diabetes$sex <- factor(diabetes$sex, levels = c(1, 2), labels = c("Male", "Female"))

cat("\nPrepared Dataset Summary:\n")
summary(diabetes)

# ------------------------------------------------------------------------------
# 2. Question 2: Explore the Data
# ------------------------------------------------------------------------------

cat("\n=== Exploratory Data Analysis ===\n")

# Descriptive statistics for outcome and key continuous predictors
cat("\nOutcome (Progression) Descriptive Stats:\n")
print(describe_numeric(diabetes$progression))

cat("\nBMI Descriptive Stats:\n")
print(describe_numeric(diabetes$bmi))

cat("\nAge Descriptive Stats:\n")
print(describe_numeric(diabetes$age))

cat("\nBlood Pressure Descriptive Stats:\n")
print(describe_numeric(diabetes$bp))

# Plot distributions using Base R
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

# Histogram of Disease Progression
hist(diabetes$progression,
     main = "Distribution of Progression",
     xlab = "Progression Score",
     col = "steelblue",
     border = "white")

# Histogram of BMI
hist(diabetes$bmi,
     main = "Distribution of BMI",
     xlab = "Body Mass Index (BMI)",
     col = "darkorange",
     border = "white")

# Boxplot of Progression
boxplot(diabetes$progression,
        main = "Boxplot of Progression",
        ylab = "Progression Score",
        col = "lightblue")

# Boxplot of Progression by Sex
boxplot(progression ~ sex,
        data = diabetes,
        main = "Progression by Sex",
        xlab = "Sex",
        ylab = "Progression Score",
        col = c("lightcyan", "pink"))

par(mfrow = c(1, 1))

# ------------------------------------------------------------------------------
# 3. Question 3: Investigate Relationships
# ------------------------------------------------------------------------------

cat("\n=== Correlation Analysis ===\n")

# Compute Pearson Correlation Matrix across continuous variables
numeric_vars <- c("age", "bmi", "bp", "s1", "s2", "s3", "s4", "s5", "s6", "progression")
cor_matrix <- cor(diabetes[, numeric_vars])

cat("\nCorrelation Matrix (Rounded to 3 decimal places):\n")
print(round(cor_matrix, 3))

# Bivariate Scatterplots with Linear Regression Trends
par(mfrow = c(1, 2))

# Progression vs BMI
plot(diabetes$bmi, diabetes$progression,
     main = "Progression vs. BMI",
     xlab = "BMI",
     ylab = "Diabetes Progression",
     pch = 19,
     col = rgb(0.2, 0.4, 0.8, 0.5))
abline(lm(progression ~ bmi, data = diabetes), col = "red", lwd = 2)

# Progression vs Blood Pressure
plot(diabetes$bp, diabetes$progression,
     main = "Progression vs. BP",
     xlab = "Blood Pressure (BP)",
     ylab = "Diabetes Progression",
     pch = 19,
     col = rgb(0.8, 0.3, 0.2, 0.5))
abline(lm(progression ~ bp, data = diabetes), col = "darkred", lwd = 2)

par(mfrow = c(1, 1))

# ------------------------------------------------------------------------------
# 4. Question 4: Build and Interpret Regression Models
# ------------------------------------------------------------------------------

cat("\n=== Linear Regression Models ===\n")

# Model 1: Simple Linear Regression (BMI)
m1_simple <- lm(progression ~ bmi, data = diabetes)
cat("\n--- Model 1: Simple Linear Regression (bmi) ---\n")
print(summary(m1_simple))

# Model 2: Multiple Linear Regression (BMI + BP + S5)
m2_multi <- lm(progression ~ bmi + bp + s5, data = diabetes)
cat("\n--- Model 2: Parsimonious Multiple Regression (bmi + bp + s5) ---\n")
print(summary(m2_multi))

# Model 3: Full Linear Regression (All Predictors)
m3_full <- lm(progression ~ age + sex + bmi + bp + s1 + s2 + s3 + s4 + s5 + s6, data = diabetes)
cat("\n--- Model 3: Full Multiple Regression (All Variables) ---\n")
print(summary(m3_full))

# Model Comparison Summary
cat("\nModel Fit Comparison (Adjusted R-squared):\n")
cat("Model 1 (BMI only):      ", round(summary(m1_simple)$adj.r.squared, 4), "\n")
cat("Model 2 (BMI + BP + S5): ", round(summary(m2_multi)$adj.r.squared, 4), "\n")
cat("Model 3 (Full Model):    ", round(summary(m3_full)$adj.r.squared, 4), "\n")

# ------------------------------------------------------------------------------
# 5. Question 5: Evaluate the Analysis
# ------------------------------------------------------------------------------

cat("\n=== Model Diagnostics & Evaluation ===\n")

# Confidence Intervals for Model 2 Coefficients
cat("\n95% Confidence Intervals for Model 2:\n")
print(confint(m2_multi))

# Multicollinearity Analysis using VIF on Full Model
cat("\nVariance Inflation Factors (VIF) for Full Model:\n")
vif_results <- manual_vif(m3_full)
print(round(vif_results, 2))

# Residual Diagnostic Plots for Model 2
par(mfrow = c(2, 2))
plot(m2_multi)
par(mfrow = c(1, 1))

cat("\n=== R Analysis Pipeline Complete ===\n")