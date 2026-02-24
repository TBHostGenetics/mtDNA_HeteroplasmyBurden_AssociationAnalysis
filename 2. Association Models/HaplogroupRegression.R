library(MASS)
library(dplyr)
library(forestplot)

# Make sure categorical variables are factors
analysis_data$sex <- as.factor(analysis_data$sex)
analysis_data$ancestry <- as.factor(analysis_data$ancestry)
analysis_data$macrohaplogroup <- as.factor(analysis_data$macrohaplogroup)

# Set reference level for haplogroup
analysis_data$macrohaplogroup <- relevel(analysis_data$macrohaplogroup, ref = "L0d1")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 1: Association between total number of heteroplasmies and haplogroup
# ══════════════════════════════════════════════════════════════════════════════

model_total <- glm.nb(total_heteroplasmy_5_95 ~ age + sex + macrohaplogroup + log(mean_read_depth) + cohort,
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
# STEP 2: Association between number of low level heteroplasmies and haplogroup
# ══════════════════════════════════════════════════════════════════════════════

model_low <- glm.nb(low_heteroplasmy_5_10 ~ age + sex + macrohaplogroup + log(mean_read_depth) + cohort,
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

# ═══════════════════════════════════════════════════════════════════════════════════════
# STEP 3: Association between number of intermediate level heteroplasmies and haplogroup
# ═══════════════════════════════════════════════════════════════════════════════════════

model_int <- glm.nb(intermediate_heteroplasmy_10_95 ~ age + sex + macrohaplogroup + log(mean_read_depth) + cohort,
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

# ══════════════════════════════════════════════════════════════════
# STEP 4: Association between number of homoplasmies and haplogroup
# ══════════════════════════════════════════════════════════════════

model_homo <- glm.nb(homoplasmy_gt95 ~ age + sex + macrohaplogroup + log(mean_read_depth) + cohort,
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

# Total number of heteroplasmic variants - results 
# Create data frame
df <- data.frame(
  term = c("Haplogroup L0a1", "Haplogroup L0a2", "Haplogroup L0d2",
           "Haplogroup L0d3", "Haplogroup L0f1", "Haplogroup L1",
           "Haplogroup L2", "Haplogroup L3", "Haplogroup L4",
           "Haplogroup L5"),
  IRR = c(5.35, 0.52, 1.76, 0.93, 4.22, 1.07, 2.71, 1.86, 1.34, 0.81),
  lower = c(3.54, 0.17, 1.29, 0.51, 1.78, 0.53, 1.86, 1.30, 0.51, 0.33),
  upper = c(8.17, 1.30, 2.40, 1.65, 10.56, 2.08, 3.96, 2.66, 3.35, 1.80),
  padj = c(4.44E-14, 0.42, 1.99E-03, 0.84, 3.83E-03, 0.84, 1.67E-06, 2.55E-03, 0.73, 0.78)
)

haplo_order <- c("L0a1","L0a2","L0d2","L0d3","L0f1","L1","L2","L3","L4","L5")

df <- df %>%
  mutate(label = gsub("Haplogroup ", "", term),
         label = factor(label, levels = haplo_order)) %>%
  arrange(label)

tabletext <- cbind(
  c("Haplogroup", as.character(df$label)),
  c("IRR (95% CI)",
    sprintf("%.2f (%.2f–%.2f)", df$IRR, df$lower, df$upper)),
  c("P-value",
    ifelse(df$padj < 0.05, 
           sprintf("%.2E", df$padj),
           sprintf("%.2f", df$padj)))
)


forestplot(
  labeltext = tabletext,
  mean  = c(NA, df$IRR),
  lower = c(NA, df$lower),
  upper = c(NA, df$upper),
  is.summary = c(TRUE, rep(FALSE, nrow(df))),
  zero = 1,
  xlog = TRUE,
  boxsize = 0.1,
  vertices = TRUE,
  fn.ci_norm = "fpDrawCircleCI",
  xticks = c(0.1, 0.5, 1, 11),
  lineheight = unit(6, "mm"),
  col = fpColors(box = "coral2", line = "coral2", zero = "grey"),
  xlab = "Incidence Rate Ratio (IRR)"
)

# Total number of homoplasmic variants - results
# Create data frame
df <- data.frame(
  term = c("Haplogroup L0a1", "Haplogroup L0a2", "Haplogroup L0d2",
           "Haplogroup L0d3", "Haplogroup L0f1", "Haplogroup L1",
           "Haplogroup L2", "Haplogroup L3", "Haplogroup L4",
           "Haplogroup L5"),
  IRR = c(0.35, 1.10, 1.02, 0.74, 2.33, 0.94, 0.58, 0.59, 1.06, 1.93),
  lower = c(0.24, 0.79, 0.89, 0.55, 1.56, 0.68, 0.46, 0.48, 0.67, 1.48),
  upper = c(0.49, 1.51, 1.18, 0.97, 3.40, 1.27, 0.73, 0.72, 1.61, 2.51),
  padj = c(2.14E-07, 0.84, 0.84, 0.08, 5.99E-05, 0.84, 1.69E-05, 2.30E-06, 0.84, 6.07E-06)
)

haplo_order <- c("L0a1","L0a2","L0d2","L0d3","L0f1","L1","L2","L3","L4","L5")

df <- df %>%
  mutate(label = gsub("Haplogroup ", "", term),
         label = factor(label, levels = haplo_order)) %>%
  arrange(label)

tabletext <- cbind(
  c("Haplogroup", as.character(df$label)),
  c("IRR (95% CI)",
    sprintf("%.2f (%.2f–%.2f)", df$IRR, df$lower, df$upper)),
  c("P-value",
    ifelse(df$padj < 0.05, 
           sprintf("%.2E", df$padj),
           sprintf("%.2f", df$padj)))
)


forestplot(
  labeltext = tabletext,
  mean  = c(NA, df$IRR),
  lower = c(NA, df$lower),
  upper = c(NA, df$upper),
  is.summary = c(TRUE, rep(FALSE, nrow(df))),
  zero = 1,
  xlog = TRUE,
  boxsize = 0.1,
  vertices = TRUE,
  fn.ci_norm = "fpDrawCircleCI",
  clip = c(0.1, 16),
  xticks = c(0.25, 0.5, 1, 2, 4),
  lineheight = unit(6, "mm"),
  col = fpColors(box = "cadetblue3", line = "cadetblue3", zero = "grey"),
  xlab = "Incidence Rate Ratio (IRR)",
)

