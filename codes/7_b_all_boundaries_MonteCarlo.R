library(ggplot2)

# Define the folder where the .rds files are stored
folder_name <- "rds_files"

# Initialize a list to store the inferred boundaries
boundary_estimates_list <- list()

# Get the list of all inferred boundary files in the folder
inferred_files <- list.files(folder_name, pattern = "inferred_boundary_OU_.*\\.rds", full.names = TRUE)

# Loop through the files and read the inferred boundary estimates
for (file in inferred_files) {
  inferred_results <- readRDS(file)
  boundary_estimates_list[[file]] <- inferred_results$boundary_est
}

# Combine all inferred boundary estimates into a matrix for easier percentile computation
boundary_matrix <- do.call(rbind, boundary_estimates_list)

# Compute the 10th and 90th percentiles across all samples
percentile_10 <- apply(boundary_matrix, 2, function(x) quantile(x, 0.10))
percentile_90 <- apply(boundary_matrix, 2, function(x) quantile(x, 0.90))

# Define the real boundary
real_boundary <- readRDS(file.path(folder_name, "boundary_OU.rds"))

# Define the timestep and calculate time points
timestep <- 1
time_points <- (0:(length(real_boundary) - 1)) * timestep

# Define the output PDF filename
pdf_filename <- "scenario_plots/all_boundaries_density.pdf"


plot(NULL, xlim = c(0, max(time_points)), ylim = range(c(percentile_10, percentile_90, real_boundary)),
     xlab = "Time", ylab = "Value")

# Add shaded area between 10th and 90th percentiles
polygon(c(time_points, rev(time_points)),
        c(percentile_10, rev(percentile_90)), col = rgb(0, 0, 1, alpha = 0.2), border = NA)

# Add the real boundary in red
lines(time_points, real_boundary, col = "red", lwd = 2)
