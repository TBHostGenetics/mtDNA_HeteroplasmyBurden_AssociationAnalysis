library(MASS)
library(dplyr)
library(forestplot)

# Make sure categorical variables are factors
analysis_data$sex <- as.factor(analysis_data$sex)
analysis_data$ancestry <- as.factor(analysis_data$ancestry)
analysis_data$macrohaplogroup <- as.factor(analysis_data$macrohaplogroup)

# Set reference level for ancestry
analysis_data$ancestry <- relevel(analysis_data$ancestry, ref = "Non-KhoeSan African")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 1: Association between total number of heteroplasmies and ancestry
# ══════════════════════════════════════════════════════════════════════════════

model_total <- glm.nb(total_heteroplasmy_5_95 ~ age + sex + ancestry + log(mean_read_depth) + cohort,
                      data = analysis_data)

# View full results
summary(model_total)

# Extract coefficients with confidence intervals
coef_summary <- summary(model_total)$coefficients

# Calculate Incidence Rate Ratios (IRR) and 95% CI
irr <- exp(coef(model_total))
irr_ci <- exp(confint(model_total))

results_total <- data.frame(
  variable = names(coef(model_total)),
  coefficient = coef(model_total),
  IRR = irr,
  CI_lower = irr_ci[,1],
  CI_upper = irr_ci[,2],
  p_value = coef_summary[,4]
)

# Apply FDR correction (Benjamini-Hochberg)
results_total$p_adjusted_FDR <- p.adjust(results_total$p_value, method = "BH")

print(results_total)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 2: Association between number of low level heteroplasmies and ancestry
# ══════════════════════════════════════════════════════════════════════════════

model_low <- glm.nb(low_heteroplasmy_5_10 ~ age + sex + ancestry + log(mean_read_depth) + cohort,
                    data = analysis_data)

# View full results
summary(model_low)

# Extract coefficients with confidence intervals
coef_summary <- summary(model_low)$coefficients

# Calculate Incidence Rate Ratios (IRR) and 95% CI
irr <- exp(coef(model_low))
irr_ci <- exp(confint(model_low))

results_low <- data.frame(
  variable = names(coef(model_low)),
  coefficient = coef(model_low),
  IRR = irr,
  CI_lower = irr_ci[,1],
  CI_upper = irr_ci[,2],
  p_value = coef_summary[,4]
)

# Apply FDR correction (Benjamini-Hochberg)
results_low$p_adjusted_FDR <- p.adjust(results_low$p_value, method = "BH")

print(results_low)

# ═════════════════════════════════════════════════════════════════════════════════════
# STEP 3: Association between number of intermediate level heteroplasmies and ancestry
# ═════════════════════════════════════════════════════════════════════════════════════

model_int <- glm.nb(intermediate_heteroplasmy_10_95 ~ age + sex + ancestry + log(mean_read_depth) + cohort,
                    data = analysis_data)

# View full results
summary(model_int)

# Extract coefficients with confidence intervals
coef_summary <- summary(model_int)$coefficients

# Calculate Incidence Rate Ratios (IRR) and 95% CI
irr <- exp(coef(model_int))
irr_ci <- exp(confint(model_int))

results_int <- data.frame(
  variable = names(coef(model_int)),
  coefficient = coef(model_int),
  IRR = irr,
  CI_lower = irr_ci[,1],
  CI_upper = irr_ci[,2],
  p_value = coef_summary[,4]
)

# Apply FDR correction (Benjamini-Hochberg)
results_int$p_adjusted_FDR <- p.adjust(results_int$p_value, method = "BH")

print(results_int)

# ════════════════════════════════════════════════════════════════
# STEP 4: Association between number of homoplasmies and ancestry
# ════════════════════════════════════════════════════════════════

model_homo <- glm.nb(homoplasmy_gt95 ~ age + sex + ancestry + log(mean_read_depth) + cohort,
                     data = analysis_data)

# View full results
summary(model_homo)

# Extract coefficients with confidence intervals
coef_summary <- summary(model_homo)$coefficients

# Calculate Incidence Rate Ratios (IRR) and 95% CI
irr <- exp(coef(model_homo))
irr_ci <- exp(confint(model_homo))

results_homo <- data.frame(
  variable = names(coef(model_homo)),
  coefficient = coef(model_homo),
  IRR = irr,
  CI_lower = irr_ci[,1],
  CI_upper = irr_ci[,2],
  p_value = coef_summary[,4]
)

# Apply FDR correction (Benjamini-Hochberg)
results_homo$p_adjusted_FDR <- p.adjust(results_homo$p_value, method = "BH")

print(results_homo)

# ════════════════════════════════════════════════════════════════
# STEP 5: Create forest plot of results
# ════════════════════════════════════════════════════════════════

# Create data frame with ancestry results from all four models
ancestry_data <- data.frame(
  model = c("Total heteroplasmic burden",
            "Number of low level heteroplasmies", 
            "Number of high level heteroplasmies",
            "Total homoplasmic burden"),
  IRR = c(0.59, 0.71, 0.55, 1.29),
  CI_lower = c(0.46, 0.50, 0.41, 1.21),
  CI_upper = c(0.76, 1.00, 0.74, 1.49),
  p_value = c(4.89E-05, 0.05, 6.28E-05, 3.95E-04), 
  p_adj = c(2.94E-04, 0.32, 3.77E-04, 2.37E-03)
)

# Add significance indicator (p < 0.05)
ancestry_data$significant <- ancestry_data$p_value < 0.05

# Create formatted text for the table
tabletext <- cbind(
  c("Outcome", ancestry_data$model),
  c("IRR", sprintf("%.2f", ancestry_data$IRR)),
  c("95% CI", sprintf("(%.2f - %.2f)", ancestry_data$CI_lower, ancestry_data$CI_upper)),
  c("P-value", ifelse(ancestry_data$p_adj < 0.05, 
                      sprintf("%.2E", ancestry_data$p_adj),
                      sprintf("%.2f", ancestry_data$p_adj)))
)

# Replace the header NA with actual text
tabletext[1, 4] <- "P-value"

# Create the forest plot
forestplot(
  tabletext,
  mean = c(NA, ancestry_data$IRR),
  lower = c(NA, ancestry_data$CI_lower),
  upper = c(NA, ancestry_data$CI_upper),
  is.summary = c(TRUE, rep(FALSE, nrow(ancestry_data))),
  clip = c(0.3, 2),
  xticks = c(0.3, 1, 2),
  xlog = FALSE,
  zero = 1,
  boxsize = 0.03,
  fn.ci_norm = "fpDrawCircleCI",
  col = fpColors(box = "darkolivegreen3", line = "darkolivegreen3", zero = "grey"),
  vertices = TRUE,
  xlab = "Incidence Rate Ratio (IRR)"
)


