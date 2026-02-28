library(ggplot2)

save_path = "rational"

# Function to compute proportions of non-inclusions from saved .rds files
compute_proportions_from_files <- function(M = 100, delta) {
  folder_name <- "rds_files"
  file_name <- paste0("boundary_", save_path, ".rds")
  file_path <- file.path(folder_name, file_name)
  
  # Initialize vectors to store results
  non_inclusions_list <- vector("list", M)
  real_boundary <- readRDS(file_path)
  cat("real_boundary =", real_boundary, "\n")
  
  # Loop through each trial index
  for (trial_index in 1:M) {
    x_sample_file <- file.path(folder_name, paste0("path_",save_path,"_", trial_index, ".rds"))
    inferred_boundary_file <- file.path(folder_name, paste0("inferred_boundary_",save_path,"_", trial_index, ".rds"))
    
    # Load the data
    X_sample <- readRDS(x_sample_file)
    inferred_results <- readRDS(inferred_boundary_file)
    if (trial_index==100){
    }
    
    # Calculate non-inclusion points
    non_inclusions <- integer(length(real_boundary))
    for (j in 1:length(real_boundary)) {
      if (real_boundary[j] < inferred_results$lower_bound[j] || real_boundary[j] > inferred_results$upper_bound[j]) {
        non_inclusions[j] <- 1
      }
    }
    
    # Store the results
    non_inclusions_list[[trial_index]] <- non_inclusions
  }
  
  # Compute proportions
  proportions <- Reduce("+", non_inclusions_list) / M
  
  # Define time points (scaled to range from 0 to 1)
  timestep <- delta
  time_points <- (0:(length(proportions) - 1)) * timestep
  time_points_scaled <- time_points / max(time_points)
  
  return(list(proportions = proportions, time_points_scaled = time_points_scaled))
}

# Set parameters
M <- 200
delta <- 1

# Compute proportions from saved files
results <- compute_proportions_from_files(M = M, delta = delta)
proportions <- results$proportions
time_points_scaled <- results$time_points_scaled

# Define alpha and calculate q_alpha
alpha <- 0.1
q_alpha <- qnorm(1 - alpha / 2)

# Calculate dotted lines values (scaled to range from 0 to 1)
dotted_lines <- alpha + c(-1, 1) * q_alpha * sqrt(alpha * (1 - alpha) / M)

# Plot proportions
plot(time_points_scaled, proportions, type = "l", col = "red", xlim = c(0, 1), ylim = c(0, 0.4),
     xlab = "Time", ylab = "Proportion of Non-inclusions")

# Add dashed line for alpha
abline(h = alpha, col = "black", lty = 2)

# Add dotted lines
abline(h = dotted_lines, col = "black", lty = 3)
