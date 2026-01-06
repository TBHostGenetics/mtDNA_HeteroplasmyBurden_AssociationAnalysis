#########################################################
# Categorize variants by variant allele frequency (VAF) #
#########################################################

# Load required libraries
library(vcfR)
library(dplyr)
library(tidyr)

# Read in VCF and covariate file
vcf <- read.vcfR("TargetedSequencingData_allSamples_PASS_filteredSamples_ExcludeBlacklistSNPs.recode.vcf")
TargetdSeq <- read.table("TargetedSeq_covariates.txt", header=TRUE)

# Extract per-sample AF from FORMAT field
af_matrix <- extract.gt(vcf, element = "AF")

# Function to handle "." and convert to numeric VAF
clean_af <- function(af_string) {
  if (is.na(af_string) || af_string == "." || af_string == "") {
    return(0)  # Reference allele
  }
  return(as.numeric(af_string))
}

# Apply to entire matrix
vaf_matrix <- apply(af_matrix, c(1, 2), clean_af)

# Function to categorize VAF
categorize_vaf <- function(vaf) {
  case_when(
    vaf == 0 ~ "reference",  # Reference allele - exclude from counts
    vaf > 0 & vaf < 0.05 ~ "very_low",  # Below your threshold
    vaf >= 0.05 & vaf < 0.10 ~ "low_heteroplasmy",
    vaf >= 0.10 & vaf < 0.95 ~ "intermediate_heteroplasmy",
    vaf >= 0.95 ~ "homoplasmy",
    TRUE ~ "other"
  )
}

# Categorize all VAF values
vaf_category_matrix <- apply(vaf_matrix, c(1, 2), categorize_vaf)

# Count variants per category per sample (excluding reference alleles)
sample_names <- colnames(vaf_matrix)
variant_counts <- data.frame(
  sample_id = sample_names,
  low_heteroplasmy_5_10 = colSums(vaf_category_matrix == "low_heteroplasmy", na.rm = TRUE),
  intermediate_heteroplasmy_10_95 = colSums(vaf_category_matrix == "intermediate_heteroplasmy", na.rm = TRUE),
  total_heteroplasmy_5_95 = colSums(vaf_category_matrix == "low_heteroplasmy", na.rm = TRUE) + 
    colSums(vaf_category_matrix == "intermediate_heteroplasmy", na.rm = TRUE),
  homoplasmy_gt95 = colSums(vaf_category_matrix == "homoplasmy", na.rm = TRUE),
  very_low_heteroplasmy_lt5 = colSums(vaf_category_matrix == "very_low", na.rm = TRUE),
  reference_alleles = colSums(vaf_category_matrix == "reference", na.rm = TRUE),
  total_variants = colSums(vaf_category_matrix != "reference", na.rm = TRUE)
)

# Create final counts for your analysis
variant_counts_final <- variant_counts %>%
  dplyr::select(sample_id, 
         low_heteroplasmy_5_10, 
         intermediate_heteroplasmy_10_95, 
         total_heteroplasmy_5_95,
         homoplasmy_gt95)

# Write to file 
write.csv(variant_counts_final, "TargetedSeq_VariantCountsByVAF.csv", row.names = FALSE)
head(variant_counts_final)

# Combine variant counts dataframe with covariate dataframe
analysis_data <- merge(TargetdSeq, variant_counts_final, by="sample_id")
