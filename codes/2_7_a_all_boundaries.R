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

# Define the real boundary
real_boundary <- readRDS(file.path(folder_name, "boundary_OU.rds"))

# Define the timestep
timestep <- 1

# Calculate time points based on the timestep and length of the boundary estimates
time_points <- (0:(length(boundary_estimates_list[[1]]) - 1)) * timestep


# Save the plot as a PDF with specified dimensions
plot(NULL, xlim = c(0, max(time_points)), ylim = range(sapply(boundary_estimates_list, range)),
     xlab = "Time", ylab = "Value")

# Add each boundary estimate to the plot
for (boundary_estimate in boundary_estimates_list) {
  lines(time_points, boundary_estimate, col = rgb(0, 0, 1, alpha = 0.1))
}

# Add the real boundary in red
lines(time_points, real_boundary, col = "red", lwd = 2)

