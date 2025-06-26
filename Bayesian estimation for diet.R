# This script implements Bayesian dietary modeling using JAGS.
# To reproduce the results presented in the supplementary text, users must manually construct Excel input files based on Supplementary Data 2 (isotopic values of food sources) and Supplementary Data 3 (isotopic values and priors for human individuals). 
# The Excel file should contain four sheets in the following order:
#   (1) human isotopic values, "sheet = "data"
#   (2) corresponding standard deviations, sheet = "error"
#   (3) Dirichlet prior parameters, sheet = "prior", and 
#   (4) isotopic values of dietary sources, sheet = "food".
#
# Note: Because isotopic values of both humans and food sources vary across time periods and regions, each group of human individuals must be matched with a specific set of dietary source values. These input files must therefore be constructed separately for each context (e.g., Anglo-Saxon England, Anglo-Saxon Scotland, Anglo-Saxon Wales). 
# The R code remains unchanged across cases; only the value of N_h (number of individuals) should be updated accordingly. 

# We provide a template dataset: "Data for running Bayesian diet coding.xlsx", which demonstrates the required data structure for running the coding. The template is based on Palaeolithic samples from Wales (6 individuals). If you want to use this dataset to run the model, please change line 51 from: N_h <- 500 to: N_h <- 6

# If you require the full set of Excel files used for running this coding, please contact us directly.

library(rjags)
library(R2jags)
library(readxl)

# Load Excel file
file_path <- "Data for running Bayesian diet coding.xlsx"

# Load observed isotopic values
obsvn_ <- as.matrix(read_excel(file_path, sheet = "data", col_names = FALSE))

# Load observation errors (standard deviations)
obsvn_error <- as.matrix(read_excel(file_path, sheet = "error", col_names = FALSE))

# Load Dirichlet prior parameters for dietary proportions（Note that alpha here corresponds to theta in the Methods section of the main text, which represents the key parameters of interest in our analysis)
a_alpha <- as.matrix(read_excel(file_path, sheet = "prior", col_names = FALSE))

# Load isotopic means and uncertainties for food sources
mI_sI_df <- as.matrix(read_excel(file_path, sheet = "food", col_names = FALSE))

# Reorganize m_I and s_I arrays into (5 food groups × 2 fractions × 2 isotopes)
m_I <- array(c(
  mI_sI_df[seq(1, 10, by=2), 1],  # 13C Protein
  mI_sI_df[seq(2, 10, by=2), 1],  # 13C Energy
  mI_sI_df[seq(1, 10, by=2), 3],  # 15N Protein
  mI_sI_df[seq(2, 10, by=2), 3]   # 15N Energy
), dim = c(5,2,2))

s_I <- array(c(
  mI_sI_df[seq(1, 10, by=2), 2],  # 13C Protein uncertainty
  mI_sI_df[seq(2, 10, by=2), 2],  # 13C Energy uncertainty
  mI_sI_df[seq(1, 10, by=2), 4],  # 15N Protein uncertainty
  mI_sI_df[seq(2, 10, by=2), 4]   # 15N Energy uncertainty
), dim = c(5,2,2))

# Number of human individuals (e.g., 500 samples，should be updated according input files)
N_h <- 500

# JAGS model specification
model_string <- "
model{

# Likelihood
  for (h_ in 1:N_h){  # individuals
    for (k_ in 1:N_k){  # isotopes
      for (j_ in 1:N_j){  # fractions (e.g., protein, energy)
        for (i_ in 1:N_i){  # food sources

          source.component.value_[h_,i_,j_,k_] <- alpha_[h_,i_] * C_[i_,j_] * (I_[i_,j_,k_] + T_[k_])
          source.component.weight_[h_,i_,j_,k_] <- alpha_[h_,i_] * C_[i_,j_]

        }
        component.contrib_[h_, j_, k_] <- W_[j_,k_] * sum(source.component.value_[h_,1:N_i,j_,k_]) / sum(source.component.weight_[h_,1:N_i,j_,k_])
      }
      H_[h_,k_] <- sum(component.contrib_[h_,1:N_j,k_]) / sum(W_[1:N_j,k_])

      # Observation model
      obsvn_[h_,k_] ~ dnorm(H_[h_,k_], tau_H[h_,k_])
      tau_H[h_,k_] <- pow(obsvn_error[h_,k_], -2)
    }
  }


  # Priors
  for (h_ in 1:N_h){
    alpha_[h_,1:N_i] ~ ddirich(a_alpha[h_,1:N_i])
  }

  # Priors for isotopic values (I_)
  I_[1,1,1] ~ dnorm(m_I[1,1,1], tau_I[1,1,1])
  tau_I[1,1,1] <- pow(s_I[1,1,1],-2)
  I_[1,2,1] ~ dnorm(m_I[1,2,1], tau_I[1,2,1])
  tau_I[1,2,1] <- pow(s_I[1,2,1],-2)
  I_[1,1,2] ~ dnorm(m_I[1,1,2], tau_I[1,1,2])
  tau_I[1,1,2] <- pow(s_I[1,1,2],-2)
  I_[1,2,2] <- 0

  I_[2,1,1] ~ dnorm(m_I[2,1,1], tau_I[2,1,1])
  tau_I[2,1,1] <- pow(s_I[2,1,1],-2)
  I_[2,2,1] ~ dnorm(m_I[2,2,1], tau_I[2,2,1])
  tau_I[2,2,1] <- pow(s_I[2,2,1],-2)
  I_[2,1,2] ~ dnorm(m_I[2,1,2], tau_I[2,1,2])
  tau_I[2,1,2] <- pow(s_I[2,1,2],-2)
  I_[2,2,2] <- 0

  I_[3,1,1] ~ dnorm(m_I[3,1,1], tau_I[3,1,1])
  tau_I[3,1,1] <- pow(s_I[3,1,1],-2)
  I_[3,2,1] ~ dnorm(m_I[3,2,1], tau_I[3,2,1])
  tau_I[3,2,1] <- pow(s_I[3,2,1],-2)
  I_[3,1,2] ~ dnorm(m_I[3,1,2], tau_I[3,1,2])
  tau_I[3,1,2] <- pow(s_I[3,1,2],-2)
  I_[3,2,2] <- 0

  I_[4,1,1] ~ dnorm(m_I[4,1,1], tau_I[4,1,1])
  tau_I[4,1,1] <- pow(s_I[4,1,1],-2)
  I_[4,2,1] ~ dnorm(m_I[4,2,1], tau_I[4,2,1])
  tau_I[4,2,1] <- pow(s_I[4,2,1],-2)
  I_[4,1,2] ~ dnorm(m_I[4,1,2], tau_I[4,1,2])
  tau_I[4,1,2] <- pow(s_I[4,1,2],-2)
  I_[4,2,2] <- 0

  I_[5,1,1] ~ dnorm(m_I[5,1,1], tau_I[5,1,1])
  tau_I[5,1,1] <- pow(s_I[5,1,1],-2)
  I_[5,2,1] ~ dnorm(m_I[5,2,1], tau_I[5,2,1])
  tau_I[5,2,1] <- pow(s_I[5,2,1],-2)
  I_[5,1,2] ~ dnorm(m_I[5,1,2], tau_I[5,1,2])
  tau_I[5,1,2] <- pow(s_I[5,1,2],-2)
  I_[5,2,2] <- 0

  # Priors for trophic offset
  for (k_ in 1:N_k){
    T_[k_] ~ dnorm(m_T[k_], tau_T[k_])
    tau_T[k_] <- pow(s_T[k_], -2)
  }

  # Priors for concentration matrix C_
  for (j_ in 1:N_j){
    for (i_ in 1:N_i){
      C_[i_,j_] ~ dnorm(m_C[i_,j_], tau_C[i_,j_]) T(0,100)
      tau_C[i_,j_] <- pow(s_C[i_,j_], -2)
    }
  }

  # Priors for weighting matrix W_
  for (k_ in 1:N_k){
    for (j_ in 1:N_j){
      W_[j_,k_] ~ dnorm(m_W[j_,k_], tau_W[j_,k_]) T(0,100)
      tau_W[j_,k_] <- pow(s_W[j_,k_], -2)
    }
  }
}
"


# Model data input
model_data <- list(
  N_h = N_h,
  N_i = 5,
  N_j = 2,
  N_k = 2,
  a_alpha = a_alpha,
  obsvn_ = obsvn_,
  obsvn_error = obsvn_error,
  m_T = c(4.8, 5.5),
  s_T = c(0.5, 0.5),
  m_W = matrix(c(74, 100, 26, 0), nrow = 2, byrow = TRUE),
  s_W = matrix(c(4, 0.001, 4, 0.001), nrow = 2, byrow = TRUE),
  m_C = matrix(c(60,40, 60,40, 10,90, 75,25, 65,35), nrow = 5, byrow = TRUE),
  s_C = matrix(c(5,5,5,5,5,5,5,5,5,5), nrow = 5, byrow = TRUE),
  m_I = m_I,
  s_I = s_I
)

# Initial values for chains
inits <- function() {
  list(alpha_ = a_alpha / 10)
}

# Run JAGS model
jags_model <- jags.model(textConnection(model_string),
                         data = model_data,
                         inits = inits,
                         n.chains=3,
                         n.adapt=5000)

update(jags_model, 5000)  # Burn-in

params <- c("alpha_")
results <- coda.samples(jags_model, variable.names=params, n.iter=5000)

# Save results to RData
save(results, file = "results.RData")

# Extract summary
results_summary <- summary(results)

# Extract statistics and quantiles for alpha_
alpha_statistics <- results_summary$statistics
alpha_statistics <- alpha_statistics[grepl("alpha_", rownames(alpha_statistics)), ]

alpha_quantiles <- results_summary$quantiles
alpha_quantiles <- alpha_quantiles[grepl("alpha_", rownames(alpha_quantiles)), ]

# Extract components
alpha_means <- alpha_statistics[, "Mean"]
alpha_sds <- alpha_statistics[, "SD"]
alpha_naive_se <- alpha_statistics[, "Naive SE"]
alpha_time_se <- alpha_statistics[, "Time-series SE"]

alpha_q_2.5 <- alpha_quantiles[, "2.5%"]
alpha_q_25  <- alpha_quantiles[, "25%"]
alpha_q_50  <- alpha_quantiles[, "50%"]
alpha_q_75  <- alpha_quantiles[, "75%"]
alpha_q_97.5 <- alpha_quantiles[, "97.5%"]

# Reshape for output (500 samples × 5 sources)
reshape_matrix <- function(vector) {
  matrix(vector, nrow = N_h, ncol = 5, byrow = FALSE)
} 

# Combine all outputs
final_df <- cbind(
  reshape_matrix(alpha_means),
  reshape_matrix(alpha_sds),
  reshape_matrix(alpha_naive_se),
  reshape_matrix(alpha_time_se),
  reshape_matrix(alpha_q_2.5),
  reshape_matrix(alpha_q_25),
  reshape_matrix(alpha_q_50),
  reshape_matrix(alpha_q_75),
  reshape_matrix(alpha_q_97.5)
)

final_df <- as.data.frame(final_df)

# Column and row labels
stat_types <- c("Mean", "SD", "NaiveSE", "TimeSE", "Q2.5", "Q25", "Q50", "Q75", "Q97.5")
colnames(final_df) <- paste0(rep(stat_types, each = 5), "_Food_", rep(1:5, times = 9))
rownames(final_df) <- paste0("Sample_", 1:N_h)

# Save CSV
write.csv(final_df, "alpha_summary_combined.csv", row.names = TRUE)