# This script is directly applied to "Supplementary Data 3.xlsx".
# To ensure consistency with the results reported in the main text, incremental enamel samples from medieval and post-medieval periods must be excluded prior to analysis.

library(readxl)
library(dplyr)
library(purrr)
library(writexl)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(patchwork)

# Set path to the Excel data file
file_path <- "Supplementary Data 3.xlsx"

# Get all sheet names from the Excel file, excluding sheet "PA"
all_sheets <- excel_sheets(file_path)
sheets_to_process <- setdiff(all_sheets, "PA")

# Define the column indices to extract: E=5, F=6, R=18 to V=22
selected_columns <- c(5, 6, 18, 19, 20, 21, 22)

# Function to read selected columns from each sheet
read_selected_columns <- function(sheet_name) {
  df <- read_excel(file_path, sheet = sheet_name, skip = 2, col_names = FALSE)
  df_selected <- df[, selected_columns]
  df_selected$sheet_name <- sheet_name
  return(df_selected)
}

# Read and combine selected data from all sheets
merged_data <- map_df(sheets_to_process, read_selected_columns)

# Rename columns for clarity
colnames(merged_data)[1:7] <- c(
  "Date_Left",
  "Date_Right",
  "Terrestrial_herbivores",
  "Terrestrial_omnivores",
  "C3_plants",
  "Marine_fish",
  "Freshwater_fish"
)

# Generate a random year within the given date ranges
set.seed(123)
merged_data <- merged_data %>%
  mutate(
    Random_Year = floor(runif(n(), min = Date_Left, max = Date_Right + 1))
  )

# Assign time bins based on the random year
merged_data <- merged_data %>%
  mutate(
    Time_Bin_Start = case_when(
      Random_Year <= -5000 ~ floor(Random_Year / 1000) * 1000,
      Random_Year > -5000 & Random_Year <= -3800 ~ -5000,
      Random_Year > -1400 & Random_Year <= -800 ~ floor(Random_Year / 200) * 200,
      TRUE ~ floor(Random_Year / 100) * 100
    ),
    Time_Bin_End = case_when(
      Random_Year <= -5000 ~ Time_Bin_Start + 1000,
      Random_Year > -5000 & Random_Year <= -3800 ~ -3800,
      Random_Year > -1400 & Random_Year <= -800 ~ Time_Bin_Start + 200,
      TRUE ~ Time_Bin_Start + 100
    ),
    Time_Bin_Label = paste0(Time_Bin_Start, "–", Time_Bin_End),
    Time_Bin_Mid = (Time_Bin_Start + Time_Bin_End) / 2
  )

# Filter out samples outside the range of interest
merged_data <- merged_data %>%
  filter(Time_Bin_Start >= -9000, Time_Bin_Start < 1900)

# Define the food categories of interest
food_cols <- c("Terrestrial_herbivores", "Terrestrial_omnivores", "C3_plants", "Marine_fish", "Freshwater_fish")

# MLE estimation function: return mean and SD if more than one value is present
estimate_mle <- function(x) {
  x <- na.omit(x)
  if (length(x) < 2) return(c(mu = NA, sigma = NA))
  mu_hat <- mean(x)
  sigma_hat <- sd(x)
  return(c(mu = mu_hat, sigma = sigma_hat))
}

# Apply MLE estimation to all food types across all time bins
mle_results <- merged_data %>%
  select(Time_Bin_Label, all_of(food_cols)) %>%
  pivot_longer(cols = all_of(food_cols), names_to = "Food", values_to = "Value") %>%
  group_by(Time_Bin_Label, Food) %>%
  summarise(
    mu = mean(Value, na.rm = TRUE),
    sigma = sd(Value, na.rm = TRUE),
    n = sum(!is.na(Value)),
    .groups = "drop"
  ) %>%
  arrange(Time_Bin_Label, Food)

# Add numeric midpoints for plotting
mle_results <- mle_results %>%
  mutate(
    Time_Bin_Start = as.numeric(sub("–.*", "", Time_Bin_Label)),
    Time_Bin_End = as.numeric(sub(".*–", "", Time_Bin_Label)),
    Time_Bin_Mid = (Time_Bin_Start + Time_Bin_End) / 2
  ) %>%
  arrange(Time_Bin_Mid)

# Define colors and labels for food groups
food_colors <- c(
  "Terrestrial_herbivores" = "#DBEDC5",
  "Terrestrial_omnivores" = "#F3DAC0",
  "C3_plants" = "#BCD1BC",
  "Marine_fish" = "#BFC4DA",
  "Freshwater_fish"  = "#C3E2EC"
)

food_labels <- c(
  "Terrestrial_herbivores" = "Terrestrial Herbivores",
  "Terrestrial_omnivores" = "Terrestrial Omnivores",
  "C3_plants"              = "C3 Plants",
  "Marine_fish"            = "Marine Resources",
  "Freshwater_fish"        = "Freshwater Resources"
)

# Smoothed plot of dietary proportions over time
p_smooth <- ggplot(mle_results, aes(x = Time_Bin_Mid, y = mu, color = Food)) +
  geom_smooth(method = "loess", se = TRUE, aes(fill = Food), linewidth = 1.2, alpha = 0.3, span = 0.4) +
  geom_vline(xintercept = c(-4300, -2500, -750, 43, 410, 1066, 1540), linetype = "dashed", color = "grey40", size = 0.3) +
  facet_wrap(~ Food, scales = "free_y", ncol = 1, labeller = as_labeller(food_labels)) +
  scale_color_manual(values = food_colors) +
  scale_fill_manual(values = food_colors) +
  labs(title = "Time-series British Diet (smoothed)", x = "Year (negative = BC)", y = "Estimated Proportion (%)") +
  scale_x_continuous(breaks = seq(-9000, 2000, 1000), minor_breaks = seq(-9000, 2000, 100), expand = c(0, 0)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_classic(base_size = 15) +
  theme(strip.background = element_blank(), strip.text = element_text(size = 13, face = "bold"), legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  annotate("text", x = -6625, y = Inf, label = "ML", size = 3, vjust = 2) +
  annotate("text", x = -3400, y = Inf, label = "NL", size = 3, vjust = 2) +
  annotate("text", x = -1625, y = Inf, label = "BA",  size = 3, vjust = 2) +
  annotate("text", x = -354,  y = Inf, label = "IA",  size = 3, vjust = 2) +
  annotate("text", x = 226.5, y = Inf, label = "RO",  size = 3, vjust = 2) +
  annotate("text", x = 738,   y = Inf, label = "AS", size = 3, vjust = 2) +
  annotate("text", x = 1303,  y = Inf, label = "MD", size = 3, vjust = 2) +
  annotate("text", x = 1720,  y = Inf, label = "PM", size = 3, vjust = 2)

# Construct shaded block data (mu ± sigma)
mle_results_rect <- mle_results %>%
  mutate(
    xmin = Time_Bin_Start,
    xmax = Time_Bin_End,
    ymin = pmax(mu - sigma, 0),
    ymax = pmin(mu + sigma, 1),
    Food_label = food_labels[Food]
  )

# Block plot: mean and confidence interval per bin
p_block <- ggplot(mle_results_rect) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Food), alpha = 0.3) +
  geom_segment(aes(x = xmin, xend = xmax, y = mu, yend = mu, color = Food), linewidth = 0.8) +
  facet_wrap(~ Food_label, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = food_colors) +
  scale_color_manual(values = food_colors) +
  scale_x_continuous(breaks = seq(-9000, 2000, 1000), minor_breaks = seq(-9000, 2000, 100), expand = c(0, 0)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  labs(title = "Time-binned Dietary Proportions", x = "Year (BC negative)", y = "Estimated Proportion (%)") +
  theme_classic(base_size = 15) +
  theme(strip.background = element_blank(), strip.text = element_text(size = 13, face = "bold"), legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

# Combine the two main plots side by side
combined_plot <- p_block | p_smooth
ggsave("Combined_Smooth_Block.png", plot = combined_plot, width = 18, height = 22.5, dpi = 300)

# Density plot of sample temporal distribution
p_density <- ggplot(merged_data, aes(x = Random_Year)) +
  geom_jitter(aes(y = 0), height = 0.00005, size = 0.8, alpha = 0.3, color = "#BE6C6D") +
  scale_x_continuous(breaks = seq(-9000, 2000, 1000), minor_breaks = seq(-9000, 2000, 100), expand = c(0, 0)) +
  labs(x = "Year (BC negative)", y = "Density") +
  theme_classic(base_size = 15) +
  theme(strip.background = element_blank(), strip.text = element_text(size = 13, face = "bold"), legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank()) +
  geom_vline(xintercept = c(-4300, -2500, -750, 43, 410, 1066, 1540), linetype = "dashed", color = "grey40", size = 0.3) +
  annotate("text", x = -6625, y = Inf, label = "ML", size = 3, vjust = 2) +
  annotate("text", x = -3400, y = Inf, label = "NL", size = 3, vjust = 2) +
  annotate("text", x = -1625, y = Inf, label = "BA",  size = 3, vjust = 2) +
  annotate("text", x = -354,  y = Inf, label = "IA",  size = 3, vjust = 2) +
  annotate("text", x = 226.5, y = Inf, label = "RO",  size = 3, vjust = 2) +
  annotate("text", x = 738,   y = Inf, label = "AS", size = 3, vjust = 2) +
  annotate("text", x = 1303,  y = Inf, label = "MD", size = 3, vjust = 2) +
  annotate("text", x = 1720,  y = Inf, label = "PM", size = 3, vjust = 2)

# Combine plots with density panel below
final_plot <- combined_plot / p_density + 
  plot_layout(heights = c(5, 1)) + 
  plot_annotation(tag_levels = "a")
ggsave("Full_Combined_with_Density_SideBySide.png", plot = final_plot, width = 18, height = 30, dpi = 600)

# Save MLE results to Excel
write_xlsx(mle_results, path = "MLE_Results.xlsx")

# Define a plotting function for selected intervals and food groups
plot_block_segment <- function(data, xmin_limit, xmax_limit, food_filter = NULL, save = FALSE, filename = NULL) {
  filtered_data <- data %>%
    filter(Time_Bin_Mid >= xmin_limit, Time_Bin_Mid <= xmax_limit)
  if (!is.null(food_filter)) {
    filtered_data <- filtered_data %>% filter(Food %in% food_filter)
  }
  p <- ggplot(filtered_data) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Food), alpha = 0.3) +
    geom_segment(aes(x = xmin, xend = xmax, y = mu, yend = mu, color = Food), linewidth = 0.8) +
    facet_wrap(~ Food_label, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = food_colors) +
    scale_color_manual(values = food_colors) +
    scale_x_continuous(limits = c(xmin_limit, xmax_limit), breaks = seq(xmin_limit, xmax_limit, by = 100), expand = c(0, 0)) +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0.05, 0.05))) +
    labs(title = paste0("Dietary Proportions (", xmin_limit, " to ", xmax_limit, ")"), x = "Year (BC negative)", y = "Estimated Proportion (%)") +
    theme_classic(base_size = 14) +
    theme(strip.background = element_blank(), strip.text = element_text(size = 13, face = "bold"), legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
  if (save && !is.null(filename)) {
    ggsave(filename, plot = p, width = 10, height = 5 + 3 * length(unique(filtered_data$Food)), dpi = 300)
  }
  return(p)
}

# Plot specific time intervals and food groups
plot_block_segment(mle_results_rect, -3600, -2500, food_filter = c("C3_plants", "Terrestrial_herbivores", "Terrestrial_omnivores"), save = TRUE, filename = "Block_-3600_to_-2500_C3.png")
plot_block_segment(mle_results_rect, -300, 400, food_filter = c("Marine_fish", "Freshwater_fish"), save = TRUE, filename = "Block_-300_to_400_Fish.png")
plot_block_segment(mle_results_rect, 400, 1500, food_filter = c("Marine_fish", "Freshwater_fish"), save = TRUE, filename = "Block_400_to_1500_Fish.png")
plot_block_segment(mle_results_rect, 1000, 1900, food_filter = c("Marine_fish"), save = TRUE, filename = "Block_1000_to_1900_Fish.png")
plot_block_segment(mle_results_rect, -2700, -1500, food_filter = c("C3_plants"), save = TRUE, filename = "Block_-2800_to_-1500_C3.png")





