# -------------------------------------------------------------------------
# Script: Allele Frequency Estimation from Ancient DNA
# Author: Zexuan Chen
# Description: This script performs bootstrap-based estimation of derived 
#              allele frequencies across sliding or fixed temporal bins, 
#              using ancient human DNA from Britain. It processes well-sampled 
#              and sparsely-sampled SNPs separately, generates trajectory plots 
#              with confidence intervals, and annotates SNPs with dietary and 
#              functional information. The allele_long_output_cleaned.xlsx 
#              and Targeted SNPs.xlsx have been provided.
# -------------------------------------------------------------------------
# Load required libraries
library(readxl)
library(dplyr)
library(purrr)
library(writexl)

# Read all sheet names from the Excel file
sheets <- excel_sheets("allele_long_output_cleaned.xlsx")

# Load each sheet into a list, with sheet names as list names
allele_long_list_cleaned <- lapply(sheets, function(sheet) {
  read_excel("allele_long_output_cleaned.xlsx", sheet = sheet)
})
names(allele_long_list_cleaned) <- sheets

# List of SNPs to be excluded (those with insufficient sample sizes)
snps_to_remove <- c(
  "6_32629764", "11_27632440", "11_61606683",  
  "11_71169547", "17_8875399", "22_31556103"
)

# Filter out excluded SNPs
allele_long_filtered <- allele_long_list_cleaned[!names(allele_long_list_cleaned) %in% snps_to_remove]

#-----------------------------------
# Define function to create sliding time windows
#-----------------------------------
set.seed(123)  # For reproducibility

generate_sliding_bins <- function(start = -4000, end = 1300, width = 1000, step = 100) {
  mids <- seq(end - width/2, start + width/2, by = -step)
  data.frame(
    Time_Bin_Mid = mids,
    Time_Bin_Start = mids - width/2,
    Time_Bin_End = mids + width/2
  )
}

# Create sliding time bins
bin_def <- generate_sliding_bins()

#-----------------------------------
# Assign samples to all overlapping time bins using randomized dates
#-----------------------------------
process_sliding_bin_multiple <- function(df, bin_df) {
  df$Date.left. <- as.numeric(df$Date.left.)
  df$Date.right. <- as.numeric(df$Date.right.)
  df <- df[!is.na(df$Date.left.) & !is.na(df$Date.right.) & df$Date.left. < df$Date.right., ]
  
  # Generate a random sampling date within each sample's date range
  df$RandomDate <- mapply(function(l, r) runif(1, l, r), df$Date.left., df$Date.right.)
  
  # For each sample, find all time bins it overlaps with (many-to-many matching)
  matched_bin_df <- do.call(rbind, lapply(1:nrow(df), function(i) {
    t <- df$RandomDate[i]
    matched <- bin_df[t >= bin_df$Time_Bin_Start & t < bin_df$Time_Bin_End, ]
    if (nrow(matched) == 0) return(NULL)
    cbind(df[rep(i, nrow(matched)), , drop = FALSE], matched)
  }))
  
  # Add human-readable label for each bin
  matched_bin_df$Time_Bin_Label <- paste0(matched_bin_df$Time_Bin_Start, "--", matched_bin_df$Time_Bin_End)
  
  return(matched_bin_df)
}

# Apply sliding window binning to all SNPs
allele_binned_multi <- lapply(allele_long_filtered, process_sliding_bin_multiple, bin_df = bin_def)

#-----------------------------------
# Prepare derived allele table
#-----------------------------------
derived_table <- read_excel("Targeted SNPs.xlsx")

# Reformat SNP name (replace colon with underscore)
derived_table$SNP <- gsub(":\\s*", "_", derived_table[[1]])
# Reformat allele representation (replace "/" with ">")
derived_table$Allele <- gsub("/", ">", derived_table$Allele)

# Keep only the SNPs that are present in the filtered binned dataset
filtered_derived_table <- derived_table[derived_table$SNP %in% names(allele_binned_multi), ]
# Create a lookup for derived allele
snp_to_derived <- setNames(filtered_derived_table$Derived, filtered_derived_table$SNP)

#-----------------------------------
# Function to bootstrap allele frequency estimation
#-----------------------------------
estimate_freq_bootstrap <- function(df, snp_name, n_boot = 1000) {
  derived_allele <- snp_to_derived[[snp_name]]
  
  if (is.null(derived_allele) || is.na(derived_allele)) {
    warning(paste("No derived allele for SNP:", snp_name))
    return(list(summary = NULL, raw = NULL))
  }
  
  df$IsDerived <- ifelse(df$Allele == derived_allele, 1, 0)
  bin_list <- split(df, df$Time_Bin_Label)
  
  summary_list <- list()
  raw_list <- list()
  
  for (bin_name in names(bin_list)) {
    bin_df <- bin_list[[bin_name]]
    if (nrow(bin_df) < 5) next  # Skip bins with too few samples
    
    boot_freqs <- replicate(n_boot, {
      sampled <- bin_df[sample(nrow(bin_df), replace = TRUE), ]
      mean(sampled$IsDerived)
    })
    
    summary_list[[bin_name]] <- data.frame(
      Time_Bin_Label = bin_name,
      Time_Bin_Start = bin_df$Time_Bin_Start[1],
      Time_Bin_End = bin_df$Time_Bin_End[1],
      Time_Bin_Mid = bin_df$Time_Bin_Mid[1],
      N = nrow(bin_df),
      Freq = mean(boot_freqs),
      CI_lower = quantile(boot_freqs, 0.025),
      CI_upper = quantile(boot_freqs, 0.975),
      SNP = snp_name
    )
    
    raw_list[[bin_name]] <- data.frame(
      SNP = snp_name,
      Time_Bin_Label = bin_name,
      Time_Bin_Mid = bin_df$Time_Bin_Mid[1],
      Replicate = 1:n_boot,
      Freq = boot_freqs
    )
  }
  
  return(list(
    summary = do.call(rbind, summary_list),
    raw = do.call(rbind, raw_list)
  ))
}

# Run bootstrap estimation for all SNPs
allele_freq_boot_list <- mapply(
  FUN = estimate_freq_bootstrap,
  df = allele_binned_multi,
  snp_name = names(allele_binned_multi),
  SIMPLIFY = FALSE
)

#-----------------------------------
# Sort output by time (midpoint) and replicate for better downstream processing
#-----------------------------------
allele_freq_boot_list_sorted <- lapply(allele_freq_boot_list, function(x) {
  list(
    summary = x$summary[order(x$summary$Time_Bin_Mid), ],
    raw     = x$raw[order(x$raw$Time_Bin_Mid, x$raw$Replicate), ]
  )
})

# Write sorted summary tables
write_xlsx(
  lapply(allele_freq_boot_list_sorted, function(x) x$summary),
  path = "allele_freq_boot_summary_sorted.xlsx"
)

# Write sorted raw bootstrap tables
write_xlsx(
  lapply(allele_freq_boot_list_sorted, function(x) x$raw),
  path = "allele_freq_boot_raw_sorted.xlsx"
)

#-----------------------------------
# Process the 6 SNPs previously excluded due to insufficient sample size
#-----------------------------------

# Extract only the excluded SNPs
allele_long_excluded <- allele_long_list_cleaned[names(allele_long_list_cleaned) %in% snps_to_remove]

#-----------------------------------
# Define fixed (non-overlapping) time bins manually
#-----------------------------------
generate_fixed_bins <- function() {
  # Manually define bin boundaries
  breaks <- c(-4000, -3000, -1000, 0, 1300)
  
  # Compute midpoints for each bin
  mids <- head(breaks, -1) + diff(breaks) / 2
  
  data.frame(
    Time_Bin_Start = head(breaks, -1),
    Time_Bin_End   = tail(breaks, -1),
    Time_Bin_Mid   = mids
  )
}

# Create the fixed bin definitions
fixed_bin_def <- generate_fixed_bins()

#-----------------------------------
# Assign each sample to its corresponding fixed bin using a randomly sampled date
#-----------------------------------
assign_fixed_bin <- function(df, bin_df) {
  df$Date.left. <- as.numeric(df$Date.left.)
  df$Date.right. <- as.numeric(df$Date.right.)
  
  # Remove invalid or missing date intervals
  df <- df[!is.na(df$Date.left.) & !is.na(df$Date.right.) & df$Date.left. < df$Date.right., ]
  
  # Generate a random date within each sample's range
  df$RandomDate <- mapply(function(l, r) runif(1, l, r), df$Date.left., df$Date.right.)
  
  # Match each sample to a fixed bin (one-to-one)
  matched <- lapply(1:nrow(df), function(i) {
    t <- df$RandomDate[i]
    matched_bin <- bin_df[t >= bin_df$Time_Bin_Start & t < bin_df$Time_Bin_End, ]
    if (nrow(matched_bin) == 0) return(NULL)
    cbind(df[i, , drop = FALSE], matched_bin)
  })
  
  result <- do.call(rbind, matched)
  result$Time_Bin_Label <- paste0(result$Time_Bin_Start, "--", result$Time_Bin_End)
  result
}

# Apply fixed binning to all 9 excluded SNPs
allele_fixed_binned <- lapply(allele_long_excluded, assign_fixed_bin, bin_df = fixed_bin_def)

#-----------------------------------
# Create a derived allele lookup table for the excluded SNPs
#-----------------------------------
filtered_derived_table_fixed <- derived_table[derived_table$SNP %in% names(allele_fixed_binned), ]
snp_to_derived_fixed <- setNames(filtered_derived_table_fixed$Derived, filtered_derived_table_fixed$SNP)

#-----------------------------------
# Define simplified bootstrap estimator (with smaller threshold)
#-----------------------------------
estimate_freq_bootstrap_simple <- function(df, snp_name, snp_to_derived_fixed, n_boot = 1000) {
  derived_allele <- snp_to_derived_fixed[[snp_name]]
  
  if (is.null(derived_allele) || is.na(derived_allele)) {
    warning(paste("No derived allele for SNP:", snp_name))
    return(list(summary = NULL, raw = NULL))
  }
  
  # Mark derived allele matches
  df$IsDerived <- ifelse(df$Allele == derived_allele, 1, 0)
  bin_list <- split(df, df$Time_Bin_Label)
  
  summary_list <- list()
  raw_list <- list()
  
  for (bin_name in names(bin_list)) {
    bin_df <- bin_list[[bin_name]]
    
    # More lenient minimum sample threshold (<3 instead of <5)
    if (nrow(bin_df) < 3) next
    
    # Bootstrap the mean derived allele frequency
    boot_freqs <- replicate(n_boot, {
      sampled <- bin_df[sample(nrow(bin_df), replace = TRUE), ]
      mean(sampled$IsDerived)
    })
    
    # Store summary
    summary_list[[bin_name]] <- data.frame(
      SNP = snp_name,
      Time_Bin_Label = bin_name,
      Time_Bin_Start = bin_df$Time_Bin_Start[1],
      Time_Bin_End = bin_df$Time_Bin_End[1],
      Time_Bin_Mid = bin_df$Time_Bin_Mid[1],
      N = nrow(bin_df),
      Freq = mean(boot_freqs),
      CI_lower = quantile(boot_freqs, 0.025),
      CI_upper = quantile(boot_freqs, 0.975)
    )
    
    # Store raw replicate results
    raw_list[[bin_name]] <- data.frame(
      SNP = snp_name,
      Time_Bin_Label = bin_name,
      Time_Bin_Mid = bin_df$Time_Bin_Mid[1],
      Replicate = 1:n_boot,
      Freq = boot_freqs
    )
  }
  
  return(list(
    summary = do.call(rbind, summary_list),
    raw = do.call(rbind, raw_list)
  ))
}

#-----------------------------------
# Run bootstrap for all 6 SNPs
#-----------------------------------
allele_fixed_boot <- mapply(
  FUN = estimate_freq_bootstrap_simple,
  df = allele_fixed_binned,
  snp_name = names(allele_fixed_binned),
  MoreArgs = list(snp_to_derived_fixed = snp_to_derived_fixed),
  SIMPLIFY = FALSE
)

#-----------------------------------
# Sort bootstrap results by time (Time_Bin_Mid and Replicate)
#-----------------------------------
allele_freq_boot_list_fixed_sorted <- lapply(allele_fixed_boot, function(x) {
  list(
    summary = x$summary[order(x$summary$Time_Bin_Mid), ],
    raw     = x$raw[order(x$raw$Time_Bin_Mid, x$raw$Replicate), ]
  )
})

# Save sorted summary results to Excel
writexl::write_xlsx(
  lapply(allele_freq_boot_list_fixed_sorted, function(x) x$summary),
  path = "allele_freq_boot_fixed_summary_sorted.xlsx"
)

# Save sorted raw replicate results to Excel
writexl::write_xlsx(
  lapply(allele_freq_boot_list_fixed_sorted, function(x) x$raw),
  path = "allele_freq_boot_fixed_raw_sorted.xlsx"
)
