#=============================#
# GAM-Based Selection Modeling Pipeline
#=============================#
# This script performs a time-indexed generalized additive model (GAM) analysis on derived allele observations from ancient DNA data, using SNP-specific temporal scaling. The allele_long_output_cleaned.xlsx and Targeted SNPs.xlsx have been provided.

#=============================#
# 1. Load Required Packages
#=============================#
library(mgcv)
library(readxl)
library(tidyverse)
library(ggrepel)
library(sf)
library(rstudioapi)
library(ggmap)
library(ggspatial)
library(ggthemes)
library(agricolae)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(ggExtra)
library(gcookbook)
library(ggstatsplot)
library(rstantools)
library(MASS)
library(fitdistrplus)
library(maxLik)
library(bayestestR)
library(forecast)
library(tidymv)
library(ggh4x)
library(ggtext)
library(writexl)
library(gridExtra)
library(gtable)
library(grid)
library(scales) 
library(patchwork)
library(writexl)

#=============================#
# 2. Load and Preprocess SNP Data
#=============================#

# Read all sheets from the Excel file into a named list
sheets <- excel_sheets("allele_long_output_cleaned.xlsx")
allele_long_output_cleaned <- lapply(sheets, function(sheet) {
  read_excel("allele_long_output_cleaned.xlsx", sheet = sheet)
})
names(allele_long_output_cleaned) <- sheets

# Load derived allele info
derived_table <- read_excel("Targeted SNPs.xlsx")
derived_table$SNP <- gsub(":\\s*", "_", derived_table[[1]])
derived_table$Allele <- gsub("/", ">", derived_table$Allele)

# Remove previously excluded low-coverage SNPs
snps_to_remove <- c("6_32629764", "11_27632440", "11_61606683",  
                    "11_71169547", "17_8875399", "22_31556103")
allele_long_output_cleaned <- allele_long_output_cleaned[!names(allele_long_output_cleaned) %in% snps_to_remove]

# Filter derived allele table accordingly
filtered_derived_table <- derived_table[derived_table$SNP %in% names(allele_long_output_cleaned), ]
derived_info <- filtered_derived_table[, c("SNP", "Derived")]

#=============================#
# 3. Label Derived Alleles (SNPs = 1) vs. Ancestral (SNPs = 0)
#=============================#

allele_long_output_cleaned_labeled <- lapply(names(allele_long_output_cleaned), function(snp_name) {
  df <- allele_long_output_cleaned[[snp_name]]
  derived_allele <- derived_info$Derived[derived_info$SNP == snp_name]
  df$SNPs <- ifelse(df$Allele == derived_allele, 1, 0)
  return(df)
})
names(allele_long_output_cleaned_labeled) <- names(allele_long_output_cleaned)

#=============================#
# 4. Sample Dates and Assign Time Index (ti)
#=============================#

set.seed(42)  # For reproducibility

# Assign a random date from the sample's date range
allele_long_output_cleaned_labeled <- lapply(allele_long_output_cleaned_labeled, function(df) {
  if (!("Date.left." %in% names(df)) || !("Date.right." %in% names(df))) {
    stop("Missing Date.left. or Date.right.")
  }
  df$Date.left. <- as.numeric(df$Date.left.)
  df$Date.right. <- as.numeric(df$Date.right.)
  df$RandomDate <- mapply(function(left, right) runif(1, min = left, max = right), df$Date.left., df$Date.right.)
  return(df)
})

# Define SNP-specific ti start dates
add_ti_column <- function(df, snp_name) {
  ti_start <- switch(
    snp_name,
    "6_32634318" = -2200,
    "6_32261252" = -1100,
    "2_136608646" = -3600,
    -4000  # Default for others
  )
  df <- df[df$RandomDate >= ti_start, ]
  df$ti <- floor((df$RandomDate - ti_start) / 60)
  return(df)
}

# Apply ti index calculation
allele_long_output_cleaned_labeled_ti <- imap(allele_long_output_cleaned_labeled, add_ti_column)

# Merge all SNP data into one data frame
allele_long_output_cleaned_labeled_all <- bind_rows(allele_long_output_cleaned_labeled_ti)

#=============================#
# 5. Fit GAM Models for Each SNP
#=============================#

# Store results
gam_results <- list()

for (snp_name in names(allele_long_output_cleaned_labeled_ti)) {
  df <- allele_long_output_cleaned_labeled_ti[[snp_name]]
  
  # Skip if required columns are missing or invalid
  if (!("SNPs" %in% colnames(df)) || !("ti" %in% colnames(df))) next
  df <- df %>% filter(!is.na(SNPs), !is.na(ti))
  if (nrow(df) < 10 || length(unique(df$SNPs)) < 2) next
  
  # Fit GAM with linear + nonlinear components
  model <- tryCatch(
    gam(
      SNPs ~ ti + s(ti, k = 7, bs = "tp", m = 1),
      select = TRUE,
      method = "REML",
      family = binomial,
      data = df,
      control = gam.control(maxit = 1000)
    ),
    error = function(e) NULL
  )
  
  # Extract results only if successful and model includes both terms
  if (!is.null(model)) {
    coef_summary <- summary(model)$p.table
    smooth_summary <- summary(model)$s.table
    
    if ("ti" %in% rownames(coef_summary) && "s(ti)" %in% rownames(smooth_summary)) {
      gam_results[[snp_name]] <- tibble(
        SNP = snp_name,
        Estimate = coef_summary["ti", "Estimate"],
        StdError = coef_summary["ti", "Std. Error"],
        z_value = coef_summary["ti", "z value"],
        p_value = coef_summary["ti", "Pr(>|z|)"],
        edf_smooth = smooth_summary["s(ti)", "edf"],
        Chisq_smooth = smooth_summary["s(ti)", "Chi.sq"],
        pval_smooth = smooth_summary["s(ti)", "p-value"]
      )
    }
  }
}

#=============================#
# 6. Merge Results with RSID Info
#=============================#

# Combine all result rows into a single data frame
gam_summary_df <- bind_rows(gam_results)

# Merge with RSID info from the derived allele table
gam_summary_df <- gam_summary_df %>%
  left_join(filtered_derived_table[, c("SNP", "RSID")], by = "SNP")
  
#=============================#
# 7. Predict Derived Allele Trajectories from GAM Models
#=============================#

# This section fits a smoothed GAM curve (nonparametric only) for each SNP and predicts the expected allele frequency trajectory from ~ -4000 to ~800 BP using the model: SNPs ~ s(ti)

# Create a list to store predictions for each SNP
gam_predictions_list <- map2(
  allele_long_output_cleaned_labeled_ti,
  names(allele_long_output_cleaned_labeled_ti),
  function(df, snp_name) {
    
    # Check that necessary columns exist
    if (!("SNPs" %in% names(df)) || !("ti" %in% names(df))) return(NULL)
    df <- df %>% filter(!is.na(SNPs), !is.na(ti))
    
    # Skip if sample size too small or only one allele category
    if (nrow(df) < 10 || length(unique(df$SNPs)) < 2) return(NULL)
    
    # Fit a GAM model with a smooth function over time (ti)
    model <- tryCatch(
      gam(SNPs ~ s(ti, k = 7, bs = "tp"), 
          family = binomial, method = "REML", data = df),
      error = function(e) NULL
    )
    if (is.null(model)) return(NULL)
    
    # Create a ti sequence for prediction: from 0 to approx. 800 years (~ti = 80)
    ti_seq <- seq(0, floor((800 + 4000) / 60), by = 1)
    new_data <- data.frame(ti = ti_seq)
    
    # Predict logit(frequency) with standard errors
    pred <- predict(model, newdata = new_data, type = "link", se.fit = TRUE)
    
    # Transform predictions from logit to probability scale with 95% CI
    new_data <- new_data %>%
      mutate(
        predicted_prob = plogis(pred$fit),
        lower_ci = plogis(pred$fit - 1.96 * pred$se.fit),
        upper_ci = plogis(pred$fit + 1.96 * pred$se.fit),
        SNP = snp_name,
        Time = -4000 + 60 * ti  # Convert ti back to actual calendar time
      )
    
    return(new_data)
  }
)

# Combine all prediction data frames into one
gam_predictions_df <- bind_rows(gam_predictions_list)

#=============================#
# Time-Binned GAM Analysis per SNP
#=============================#
# This script performs generalized additive modeling (GAM) within each time bin for each SNP, capturing local selection signals. It builds on pre-binned data with derived allele labels and calculates both parametric and smooth estimates over time. The allele_binned_multi.xlsx can be generated by Allele Frequency Trajectories.R, and is thus not provided directly.

#=============================#
# 1. Load Binned SNP Data
#=============================#
sheet_names <- excel_sheets("allele_binned_multi.xlsx") #allele_binned_multi can be generated by Allele Frequency Trajectories.R

allele_binned_multi <- lapply(sheet_names, function(sheet) {
  read_excel("allele_binned_multi.xlsx", sheet = sheet)
})
names(allele_binned_multi) <- sheet_names

#=============================#
# 2. Label Derived Alleles (SNPs = 1) and Compute ti
#=============================#

allele_binned_multi_labeled <- lapply(names(allele_binned_multi), function(snp_name) {
  df <- allele_binned_multi[[snp_name]]
  derived_allele <- derived_info$Derived[derived_info$SNP == snp_name]
  df$SNPs <- ifelse(df$Allele == derived_allele, 1, 0)
  return(df)
})
names(allele_binned_multi_labeled) <- names(allele_binned_multi)

# Add time index (ti), adjusted by SNP-specific origin
allele_binned_multi_labeled_ti <- lapply(names(allele_binned_multi_labeled), function(snp_name) {
  df <- allele_binned_multi_labeled[[snp_name]]
  
  ti_start <- switch(
    snp_name,
    "6_32634318" = -2200,
    "6_32261252" = -1100,
    "2_136608646" = -3600,
    -4000
  )
  
  df <- df[df$RandomDate >= ti_start, ]
  
  # Use 12-year unit to allow finer time resolution
  df$ti <- floor((df$RandomDate - (-4000)) / 12)
  
  return(df)
})
names(allele_binned_multi_labeled_ti) <- names(allele_binned_multi_labeled)

# Combine into a single data frame
allele_binned_multi_labeled_ti_all <- bind_rows(allele_binned_multi_labeled_ti)

#=============================#
# 3. Fit Time-Bin Local GAM Models
#=============================#

gam_bin_results <- list()

for (snp_name in names(allele_binned_multi_labeled_ti)) {
  df_all <- allele_binned_multi_labeled_ti[[snp_name]]
  
  # Split by time bins
  df_split <- split(df_all, df_all$Time_Bin_Mid)
  
  # Fit one GAM per time bin
  for (bin_label in names(df_split)) {
    df <- df_split[[bin_label]]
    
    model <- tryCatch(
      gam(
        SNPs ~ ti + s(ti, k = 7, bs = "tp", m = 1),
        family = binomial,
        method = "ML",
        select = TRUE,
        data = df,
        control = gam.control(maxit = 1000, epsilon = 1e-8)
      ),
      error = function(e) NULL
    )
    
    if (!is.null(model)) {
      coef_summary <- summary(model)$p.table
      smooth_summary <- summary(model)$s.table
      
      if ("ti" %in% rownames(coef_summary) && "s(ti)" %in% rownames(smooth_summary)) {
        gam_bin_results[[paste0(snp_name, "_", bin_label)]] <- tibble(
          SNP = snp_name,
          Time_Bin_Mid = as.numeric(bin_label),
          Estimate = coef_summary["ti", "Estimate"],
          StdError = coef_summary["ti", "Std. Error"],
          z_value = coef_summary["ti", "z value"],
          p_value = coef_summary["ti", "Pr(>|z|)"],
          edf_smooth = smooth_summary["s(ti)", "edf"],
          Chisq_smooth = smooth_summary["s(ti)", "Chi.sq"],
          pval_smooth = smooth_summary["s(ti)", "p-value"]
        )
      }
    }
  }
}

# Combine results across all bins and SNPs
gam_bin_summary_df <- bind_rows(gam_bin_results)

# Apply SNP-specific time truncation to ensure meaningful ti starts
gam_bin_summary_df <- gam_bin_summary_df %>%
  filter(
    !(SNP == "6_32634318" & Time_Bin_Mid < -900),
    !(SNP == "6_32261252" & Time_Bin_Mid < -400),
    !(SNP == "2_136608646" & Time_Bin_Mid < -2700)
  )


#=============================#
# Exporting Data to Excel
#=============================#
output_list <- list(
  "GAM_mean" = gam_summary_df,
  "GAM_timebin" = gam_bin_summary_df
)

write_xlsx(output_list, path = "gam_model_summary_output.xlsx")