################################################################################
# Title: Convergent Cross Mapping (CCM) Analysis of Diet and Selection Coefficients
# Author: Zexuan Chen
# Purpose:
#   This script performs Convergent Cross Mapping (CCM) analyses to assess
#   the potential causal relationship between dietary patterns and time-varying
#   selection coefficients inferred from ancient DNA data.
#
#   It includes:
#     - CCM based on point-estimated diet and selection coefficients;
#     - 100 simulations of both dietary input and selection coefficients;
#     - Overlay of real and simulated CCM curves per SNP.
################################################################################

# Load required packages
library(readxl)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)

# === Load input data ===

# Read pottery residue summary (Dairy data)
sherd_summary <- read_excel("Supplymentary Data 9.xlsx", sheet = "Diet")

# Read isotopic dietary reconstruction (mean ± sd per time bin)
plot_data <- read_excel("Supplymentary Data 9.xlsx", sheet = "Dairy")

# Read selection coefficient estimates from GAM model
GAM_timebin <- read_excel("Supplymentary Data 8.xlsx", sheet = "GAM_timebin")

# Filter out two irrelevant SNPs
GAM_timebin_filtered <- GAM_timebin %>%
  filter(!RSID %in% c("rs3891176", "rs7775397")) %>%
  filter(!(RSID == "rs4988235" & Time_Bin_Mid > 100))  # Limit rs4988235 to ancient time bins

# Filter sherd data to time bins ≥ -2700
sherd_summary <- sherd_summary %>%
  filter(Time_Bin_Mid >= -2700)

# Define SNP-to-food group mapping
snp_food_map <- list(
  Dairy = c("rs4988235"),
  Marine_fish = c("rs7944926", "rs174546", "rs174594", "rs174570"),
  C3_plants = c("rs653178", "rs272872", "rs7517", "rs4073089", 
                "rs6265", "rs12401678", "rs2010410")
)

# Pivot isotopic dietary means by food group
diet_series_list <- plot_data %>%
  select(Time_Bin_Mid, Food, mean) %>%
  pivot_wider(names_from = Food, values_from = mean)

# Extract dairy proportion (from pottery)
dairy_series <- sherd_summary %>%
  select(Time_Bin_Mid, Dairy = Fraction_below)

# Merge all food group series into one table
diet_all <- full_join(diet_series_list, dairy_series, by = "Time_Bin_Mid")

# Get list of all SNPs from filtered GAM output
all_snp <- unique(GAM_timebin_filtered$RSID)

# Build time series for each SNP: diet + selection coefficient
ccm_input_list <- map(all_snp, function(snp) {
  snp_series <- GAM_timebin_filtered %>%
    filter(RSID == snp) %>%
    select(Time_Bin_Mid, Selection = Estimate)
  
  diet_group <- names(snp_food_map)[map_lgl(snp_food_map, ~ snp %in% .x)]
  if (length(diet_group) == 0) return(NULL)
  
  input <- diet_all %>%
    select(Time_Bin_Mid, Diet = !!sym(diet_group)) %>%
    inner_join(snp_series, by = "Time_Bin_Mid") %>%
    arrange(Time_Bin_Mid) %>%
    mutate(RSID = snp, Food_Group = diet_group)
  
  return(input)
})

# Name and clean the list
names(ccm_input_list) <- all_snp
ccm_input_list <- compact(ccm_input_list)

# === Real data CCM analysis ===

E <- 3  # Embedding dimension
tau <- 1  # Time delay
k <- 4  # Number of nearest neighbors
start_generation <- 5  # Start from minimum library length

# Create time-delay embedding matrix
embed_space <- function(series, E, tau) {
  N <- length(series)
  max_index <- N - (E - 1) * tau
  if (max_index <= 0) stop("Insufficient data to construct the embedding matrix")
  
  M <- matrix(NA, nrow = max_index, ncol = E)
  for (i in 1:E) {
    indices <- (1:max_index) + (E - i) * tau
    M[, i] <- series[indices]
  }
  return(M)
}

# Predict X from embedding of Y
predict_X <- function(M_Y, X, k) {
  n <- nrow(M_Y)
  predicted <- numeric(n)
  for (t in 1:n) {
    distances <- apply(M_Y, 1, function(row) sqrt(sum((row - M_Y[t, ])^2)))
    distances[t] <- Inf  # exclude self
    nearest_indices <- order(distances)[1:k]
    weights <- exp(-distances[nearest_indices])
    weights <- weights / sum(weights)
    predicted[t] <- sum(weights * X[nearest_indices])
  }
  return(predicted)
}

# Compute correlation between predicted and actual X
calculate_correlation <- function(M_Y, X, L, k) {
  if (L <= E) return(c(correlation = NA, lower = NA, upper = NA))
  pred <- predict_X(M_Y[1:L, ], X[1:L], k)
  cor_test <- cor.test(pred, X[1:L], conf.level = 0.95)
  return(c(correlation = cor_test$estimate, lower = cor_test$conf.int[1], upper = cor_test$conf.int[2]))
}

# Apply CCM for each SNP in the real dataset
ccm_results_realdata <- purrr::imap_dfr(ccm_input_list, function(df, snp) {
  df <- df %>% arrange(Time_Bin_Mid)
  Y <- df$Selection
  X <- df$Diet
  M_Y <- embed_space(Y, E, tau)
  if (is.null(M_Y)) return(NULL)
  X <- tail(X, nrow(M_Y))
  
  L_values <- start_generation:nrow(M_Y)
  res_df <- data.frame(L = L_values, Correlation = NA, Lower = NA, Upper = NA)
  for (i in seq_along(L_values)) {
    L <- L_values[i]
    res_df[i, 2:4] <- calculate_correlation(M_Y, X, L, k)
  }
  
  res_df$RSID <- snp
  res_df$Food_Group <- unique(df$Food_Group)
  return(res_df)
})

# === Plotting labels ===

# Prepare SNP labels (e.g., "rs4988235 (C>T)")
derived_table <- read_excel("D:\\杜伦大学\\博三科研\\nature 投稿相关\\Targeted SNPs.xlsx", sheet = "Sheet3")
derived_table$SNP <- gsub(":\\s*", "_", derived_table[[1]])
derived_table$Allele <- gsub("/", ">", derived_table$Allele)
facet_labels <- setNames(
  paste0(derived_table$RSID, " (", derived_table$Allele, ")"),
  derived_table$RSID
)

# === Plot real data CCM results ===
ggplot(ccm_results_realdata, aes(x = L, y = Correlation, color = Food_Group, fill = Food_Group)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, color = NA) +
  geom_line(size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.4) +
  facet_wrap(~ RSID, scales = "free_y", labeller = labeller(RSID = facet_labels)) +
  scale_color_manual(values = c("C3_plants" = "#BCD1BC", "Dairy" = "#F3DAC0", "Marine_fish" = "#BFC4DA")) +
  scale_fill_manual(values = c("C3_plants" = "#BCD1BC", "Dairy" = "#F3DAC0", "Marine_fish" = "#BFC4DA")) +
  labs(title = "CCM based on point-estimated diet and selection coefficients",
       x = "Library length (L)", y = expression(rho)) +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0)) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 8, face = "bold"),
        legend.position = "none",
        legend.title = element_blank())

# === Simulated selection coefficients ===

simulate_all_selection <- function(df, n_sim = 100) {
  df %>%
    mutate(StdError = ifelse(is.na(StdError), 0.01, StdError)) %>%
    mutate(Sim_values = map2(Estimate, StdError, ~ rnorm(n_sim, .x, .y))) %>%
    unnest_longer(Sim_values) %>%
    group_by(RSID, Time_Bin_Mid) %>%
    arrange(desc(Sim_values), .by_group = TRUE) %>%
    mutate(Sim_ID = row_number()) %>%
    ungroup() %>%
    rename(Sim_s = Sim_values) %>%
    select(RSID, Time_Bin_Mid, Sim_ID, Sim_s)
}
sim_s_all <- simulate_all_selection(GAM_timebin_filtered)

# === Simulated diet data ===

# Non-dairy food simulation
simulate_plot_diet <- function(df, n_sim = 100) {
  df %>%
    mutate(Sim_values = map2(mean, sd, ~ rnorm(n_sim, .x, .y))) %>%
    unnest_longer(Sim_values) %>%
    group_by(Food, Time_Bin_Mid) %>%
    arrange(desc(Sim_values), .by_group = TRUE) %>%
    mutate(Sim_ID = row_number()) %>%
    ungroup() %>%
    rename(Sim_diet = Sim_values) %>%
    select(Food, Time_Bin_Mid, Sim_ID, Sim_diet)
}
sim_diet_plot <- simulate_plot_diet(plot_data)

# Dairy simulation using normal approximation from binomial CI
simulate_dairy <- function(df, n_sim = 100) {
  df %>%
    mutate(sd_approx = (upper - lower) / (2 * 1.96)) %>%
    mutate(Sim_values = map2(Fraction_below, sd_approx, ~ rnorm(n_sim, .x, .y))) %>%
    unnest_longer(Sim_values) %>%
    group_by(Time_Bin_Mid) %>%
    arrange(desc(Sim_values), .by_group = TRUE) %>%
    mutate(Sim_ID = row_number()) %>%
    ungroup() %>%
    rename(Sim_diet = Sim_values) %>%
    mutate(Food = "Dairy") %>%
    select(Food, Time_Bin_Mid, Sim_ID, Sim_diet)
}
sim_diet_dairy <- simulate_dairy(sherd_summary)
sim_diet_all <- bind_rows(sim_diet_plot, sim_diet_dairy)

# Mapping between SNP and food group
snp_food_table <- tibble(
  Food = rep(names(snp_food_map), lengths(snp_food_map)),
  RSID = unlist(snp_food_map)
)

# Core Functions for Convergent Cross Mapping (CCM)
embed_space <- function(series, E = 3, tau = 1) {
  N <- length(series)
  max_index <- N - (E - 1) * tau
  if (max_index <= 0) return(NULL)
  
  M <- matrix(NA, nrow = max_index, ncol = E)
  for (i in 1:E) {
    indices <- (1:max_index) + (E - i) * tau
    M[, i] <- series[indices]
  }
  return(M)
}

predict_X <- function(M_Y, X, k = 4) {
  n <- nrow(M_Y)
  predicted <- numeric(n)
  for (t in 1:n) {
    distances <- apply(M_Y, 1, function(row) sqrt(sum((row - M_Y[t, ])^2)))
    distances[t] <- Inf
    nearest_indices <- order(distances)[1:k]
    weights <- exp(-distances[nearest_indices])
    weights <- weights / sum(weights)
    predicted[t] <- sum(weights * X[nearest_indices])
  }
  return(predicted)
}

calculate_correlation <- function(M_Y, X, L, k = 4) {
  if (L <= 3 || L > nrow(M_Y)) return(NA)
  pred <- predict_X(M_Y[1:L, ], X[1:L], k)
  cor(pred, X[1:L], use = "complete.obs")
}

# === CCM on simulated data ===

ccm_results <- list()

for (snp in unique(sim_s_all$RSID)) {
  cat("Processing SNP:", snp, "\n")
  sim_snp <- sim_s_all %>% filter(RSID == snp)
  food_group <- snp_food_table %>% filter(RSID == snp) %>% pull(Food)
  sim_food <- sim_diet_all %>% filter(Food == food_group)
  sim_merged <- inner_join(sim_snp, sim_food, by = c("Time_Bin_Mid", "Sim_ID"))
  
  sim_result <- sim_merged %>%
    group_by(Sim_ID) %>%
    group_split() %>%
    imap_dfr(function(df, sim_id) {
      M_Y <- embed_space(df$Sim_s)
      if (is.null(M_Y)) return(NULL)
      X <- tail(df$Sim_diet, nrow(M_Y))
      L_vals <- 5:nrow(M_Y)
      data.frame(
        L = L_vals,
        Correlation = map_dbl(L_vals, ~ calculate_correlation(M_Y, X, .x)),
        RSID = snp,
        Food_Group = food_group,
        Sim_ID = sim_id
      )
    })
  ccm_results[[snp]] <- sim_result
}

# Combine results
ccm_sim_all <- bind_rows(ccm_results)

# === Plot final comparison ===

ggplot() +
  geom_line(data = ccm_sim_all, 
            aes(x = L, y = Correlation, group = interaction(RSID, Sim_ID)),
            color = "#EAF3FA", size = 0.4) +
  geom_ribbon(data = ccm_results_realdata,
              aes(x = L, ymin = Lower, ymax = Upper, fill = Food_Group),
              alpha = 0.2, color = NA) +
  geom_line(data = ccm_results_realdata,
            aes(x = L, y = Correlation, color = Food_Group),
            size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.4) +
  facet_wrap(~ RSID, scales = "free_y", labeller = labeller(RSID = facet_labels)) +
  scale_color_manual(values = c("C3_plants" = "#BCD1BC", "Dairy" = "#F3DAC0", "Marine_fish" = "#BFC4DA")) +
  scale_fill_manual(values = c("C3_plants" = "#BCD1BC", "Dairy" = "#F3DAC0", "Marine_fish" = "#BFC4DA")) +
  labs(title = "CCM Results: (Mean vs Simulated)", x = "Library length (L)", y = expression(rho)) +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0)) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 8, face = "bold"),
        legend.position = "none",
        legend.title = element_blank())

ggsave("CCM_Results_Facet.png", width = 10, height = 8, units = "in", dpi = 300)