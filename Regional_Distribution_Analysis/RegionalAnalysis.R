# ========================================
# REGIONAL DISTRIBUTION ANALYSIS
# ========================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(vcfR)
library(gridExtra)

# ========================================
# FUNCTION 1: ANNOTATE POSITIONS TO GENE REGIONS
# ========================================

annotate_gene_region <- function(position) {
  # Human mtDNA genome (rCRS reference)
  
  # D-loop (wraps around)
  if (position >= 16024 | position <= 576) return("D-loop")
  
  # rRNA genes
  if (position >= 648 & position <= 1601) return("12S rRNA")
  if (position >= 1671 & position <= 3229) return("16S rRNA")
  
  # Protein-coding genes (OXPHOS complexes)
  if (position >= 3307 & position <= 4262) return("ND1")
  if (position >= 4470 & position <= 5511) return("ND2")
  if (position >= 5904 & position <= 7445) return("COX1")
  if (position >= 7586 & position <= 8269) return("COX2")
  if (position >= 8366 & position <= 8572) return("ATP8")
  if (position >= 8527 & position <= 9207) return("ATP6")
  if (position >= 9207 & position <= 9990) return("COX3")
  if (position >= 10059 & position <= 10404) return("ND3")
  if (position >= 10470 & position <= 10766) return("ND4L")
  if (position >= 10760 & position <= 12137) return("ND4")
  if (position >= 12337 & position <= 14148) return("ND5")
  if (position >= 14149 & position <= 15887) return("ND6")
  if (position >= 14747 & position <= 15887) return("CYTB")
  
  # Everything else (tRNA/intergenic)
  return("tRNA Other")
}

# ========================================
# FUNCTION 2: CALCULATE REGIONAL COUNTS
# ========================================

calculate_regional_counts <- function(vcf_file, metadata, vaf_type = "heteroplasmy") {
  
  cat("Processing VCF file:", vcf_file, "\n")
  
  # Read VCF
  vcf <- read.vcfR(vcf_file, verbose = FALSE)
  
  # Get positions
  positions <- getPOS(vcf)
  
  # Annotate regions
  regions <- sapply(positions, annotate_gene_region)
  
  # Extract VAF matrix
  af_matrix <- extract.gt(vcf, element = "AF")
  
  # Clean AF values
  clean_af <- function(af_string) {
    if (is.na(af_string) || af_string == "." || af_string == "") {
      return(0)  # Reference allele
    }
    return(as.numeric(af_string))
  }
  
  vaf_matrix <- apply(af_matrix, c(1, 2), clean_af)
  
  # Categorize by VAF type
  if (vaf_type == "heteroplasmy") {
    # Heteroplasmies: 5-95%
    variant_matrix <- apply(vaf_matrix, c(1, 2), function(vaf) {
      !is.na(vaf) & vaf >= 0.05 & vaf < 0.95
    })
  } else if (vaf_type == "homoplasmy") {
    # Homoplasmies: >95% (excluding haplogroup-defining variants)
    variant_matrix <- apply(vaf_matrix, c(1, 2), function(vaf) {
      !is.na(vaf) & vaf >= 0.95
    })
  } else if (vaf_type == "intermediate") {
    # Intermediate only: 10-95%
    variant_matrix <- apply(vaf_matrix, c(1, 2), function(vaf) {
      !is.na(vaf) & vaf >= 0.10 & vaf < 0.95
    })
  } else if (vaf_type == "low") {
    # Low-level only: 5-10%
    variant_matrix <- apply(vaf_matrix, c(1, 2), function(vaf) {
      !is.na(vaf) & vaf >= 0.05 & vaf < 0.10
    })
  }
  
  # Get sample IDs
  samples <- colnames(vaf_matrix)
  
  # Count by region for each sample
  results_list <- list()
  
  for (sample in samples) {
    
    # Create dataframe for this sample
    sample_df <- data.frame(
      position = positions,
      region = regions,
      has_variant = variant_matrix[, sample],
      stringsAsFactors = FALSE
    )
    
    # Count by region
    region_counts <- sample_df %>%
      group_by(region) %>%
      summarize(
        n_variants = sum(has_variant),
        .groups = 'drop'
      ) %>%
      mutate(sample_id = sample)
    
    results_list[[sample]] <- region_counts
  }
  
  # Combine all samples
  all_counts <- bind_rows(results_list)
  
  # Merge with metadata
  all_counts <- all_counts %>%
    left_join(metadata, by = "sample_id")
  
  return(all_counts)
}

# ========================================
# FUNCTION 3: STATISTICAL TESTING BY REGION
# ========================================

test_regional_differences <- function(regional_data, group_var, 
                                      group1_name, group2_name = NULL) {
  
  # If group2_name not specified, compare group1 vs all others
  if (is.null(group2_name)) {
    test_data <- regional_data
  } else {
    test_data <- regional_data %>%
      filter(!!sym(group_var) %in% c(group1_name, group2_name))
  }
  
  # Get unique regions
  regions <- unique(test_data$region)
  
  # Test each region
  results_list <- list()
  
  for (reg in regions) {
    
    region_data <- test_data %>% filter(region == reg)
    
    # Extract data for each group
    group1_data <- region_data %>% 
      filter(!!sym(group_var) == group1_name) %>% 
      pull(n_variants)
    
    if (!is.null(group2_name)) {
      group2_data <- region_data %>% 
        filter(!!sym(group_var) == group2_name) %>% 
        pull(n_variants)
    } else {
      group2_data <- region_data %>% 
        filter(!!sym(group_var) != group1_name) %>% 
        pull(n_variants)
    }
    
    # Calculate summary statistics
    n_group1 <- length(group1_data)
    n_group2 <- length(group2_data)
    mean_group1 <- mean(group1_data, na.rm = TRUE)
    mean_group2 <- mean(group2_data, na.rm = TRUE)
    sd_group1 <- sd(group1_data, na.rm = TRUE)
    sd_group2 <- sd(group2_data, na.rm = TRUE)
    difference <- mean_group1 - mean_group2
    
    # Perform Welch's t-test
    t_test <- tryCatch({
      t.test(group1_data, group2_data, var.equal = FALSE)
    }, error = function(e) {
      NULL
    })
    
    # Extract test statistics
    if (!is.null(t_test)) {
      t_statistic <- as.numeric(t_test$statistic)
      p_value <- t_test$p.value
      df <- as.numeric(t_test$parameter)
    } else {
      t_statistic <- NA
      p_value <- NA
      df <- NA
    }
    
    # Store results
    results_list[[reg]] <- data.frame(
      region = reg,
      n_group1 = n_group1,
      n_group2 = n_group2,
      mean_group1 = mean_group1,
      sd_group1 = sd_group1,
      mean_group2 = mean_group2,
      sd_group2 = sd_group2,
      difference = difference,
      t_statistic = t_statistic,
      df = df,
      p_value = p_value,
      stringsAsFactors = FALSE
    )
  }
  
  # Combine results
  results <- bind_rows(results_list)
  
  # Apply FDR correction
  results <- results %>%
    mutate(
      p_adjusted_FDR = p.adjust(p_value, method = "BH"),
      significant_FDR = p_adjusted_FDR < 0.05,
      # Calculate Cohen's d (effect size)
      pooled_sd = sqrt(((n_group1 - 1) * sd_group1^2 + (n_group2 - 1) * sd_group2^2) / 
                         (n_group1 + n_group2 - 2)),
      cohens_d = difference / pooled_sd
    ) %>%
    arrange(p_value)
  
  return(results)
}

# ========================================
# FUNCTION 4: CREATE SUMMARY PLOTS
# ========================================

plot_regional_barplot <- function(regional_data, group_var, 
                                  comparison_name, test_results) {
  
  # Calculate summary statistics
  summary_data <- regional_data %>%
    group_by(region, !!sym(group_var)) %>%
    summarize(
      mean_count = mean(n_variants),
      se = sd(n_variants) / sqrt(n()),
      .groups = 'drop'
    ) %>%
    mutate(
      region = factor(region, levels = c(
        "D-loop", "12S rRNA", "16S rRNA", 
        "ND1", "ND2", "COX1", "COX2", "ATP8", "ATP6", "COX3",
        "ND3", "ND4L", "ND4", "ND5", "ND6", "CYTB", "tRNA Other"
      ))
    )
  
  # Add significance labels
  sig_labels <- test_results %>%
    dplyr::select(region, p_adjusted_FDR) %>%
    mutate(
      sig_label = case_when(
        p_adjusted_FDR < 0.001 ~ "***",
        p_adjusted_FDR < 0.01 ~ "**",
        p_adjusted_FDR < 0.05 ~ "*",
        TRUE ~ ""
      )
    )
  
  summary_data <- summary_data %>%
    left_join(sig_labels, by = "region")
  
  # Get the first group value for positioning significance stars
  first_group <- unique(summary_data[[group_var]])[1]
  
  # Calculate max y value for positioning stars
  max_y <- max(summary_data$mean_count + summary_data$se, na.rm = TRUE)
  
  # Create plot
  p <- ggplot(summary_data, 
              aes(x = region, y = mean_count, fill = !!sym(group_var))) +
    
    geom_col(position = position_dodge(width = 0.8), alpha = 0.8, width = 0.75) +
    
    geom_errorbar(aes(ymin = mean_count - se, ymax = mean_count + se),
                  position = position_dodge(width = 0.8),
                  width = 0.3, size = 0.6) +
    
    # Add significance stars (only for one group to avoid duplication)
    geom_text(data = summary_data %>% 
                filter(!!sym(group_var) == first_group) %>%
                distinct(region, .keep_all = TRUE),
              aes(label = sig_label, y = mean_count + se + max_y * 0.1),
              position = position_dodge(width = 0.8),
              size = 6, color = "black", fontface = "bold") +
    
    labs(
      title = paste("Regional Heteroplasmy Distribution:", comparison_name),
      x = "Gene Region",
      y = "Mean Number of Heteroplasmic Variants",
      fill = group_var,
      caption = "*** FDR p < 0.001, ** p < 0.01, * p < 0.05\nError bars = standard error"
    ) +
    
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.major.x = element_blank()
    )
  
  return(p)
}

# ========================================
# EXAMPLE: TB CASES VS CONTROLS
# ========================================

cohort1_metadata <- read.table("TargetedSeq_covariates.txt", header=TRUE)

# Calculate regional counts for heteroplasmies
cohort1_regional_het <- calculate_regional_counts(
  vcf_file = "TargetedSequencingData_allSamples_PASS_filteredSamples_ExcludeBlacklistSNPs.recode.vcf",
  metadata = cohort1_metadata,
  vaf_type = "heteroplasmy"
)

# Test TB vs Control
cohort1_tests_het <- test_regional_differences(
  regional_data = cohort1_regional_het,
  group_var = "phenotype",
  group1_name = "Case",
  group2_name = "Control"
)

print(cohort1_tests_het)
write.csv(cohort1_tests_het, "cohort1_regional_heteroplasmy_tests.csv", row.names = FALSE)

# Plot
p_cohort1_het <- plot_regional_barplot(
  regional_data = cohort1_regional_het,
  group_var = "phenotype",
  comparison_name = "Cohort 1 (Active TB vs Controls)",
  test_results = cohort1_tests_het
)

print(p_cohort1_het)
ggsave("cohort1_regional_heteroplasmy.png", p_cohort1_het, 
       width = 12, height = 6, dpi = 300)
