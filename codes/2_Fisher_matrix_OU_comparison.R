library(ggplot2)
library(latex2exp)
library(reshape2)

store_matrices <- function(mu, sigma2, alpha, sample_size = 1000, T = 1, delta = 0.01) {
  n <- (T/delta)
  theta <- c(alpha, mu, sigma2)
  
  fisher_info_gradient_list <- vector("list", length = sample_size)
  
  for (i in 1:sample_size) {
    X <- <- simulate_gm_process(n_paths, n_steps, delta, x0, slope_fun_t, trend_fun_t, vol_fun_t, alpha_true, beta_true, sigma_true, subdivision)$path
    
    result <- est_OU(X, delta)
    
    fisher_info_second_derivatives_list[[i]] <- fisher_info_second_derivatives(result$alpha, result$mu, result$sigma2, X, delta,
                                                                               d2_logL_dalpha2, d2_logL_dmu2, d2_logL_dsigma2,
                                                                               d2_logL_dalphadmu, d2_logL_dalphadsigma, d2_logL_dmudsigma)
    
    fisher_info_gradient_list[[i]] <- fisher_info_gradient(result$alpha, result$mu, result$sigma2, X, delta,
                                                           d_logL_dalpha, d_logL_dmu, d_logL_dsigma2)
    
    R <- chol(fisher_info_second_derivatives_list[[i]])
    
    point <- sqrt(n) * (R %*% (c(result$alpha, result$mu, result$sigma2) - theta))
    data[i,] <- point
  }
  
  return(list(data = data,
              fisher_info_second_derivatives_list = fisher_info_second_derivatives_list,
              fisher_info_gradient_list = fisher_info_gradient_list))
}

results <- store_matrices(20, 2, 3)

# ALPHA
alpha_fisher_entry <- alpha_fisher(3, 2, 20, 26, 0.01)
i <- 1
second_derivatives <- sapply(results$fisher_info_second_derivatives_list, function(m) m[i, i])
gradients <- sapply(results$fisher_info_gradient_list, function(m) m[i, i])

df_alpha <- data.frame(value = c(second_derivatives, gradients),
                       type = rep(c("Second Derivative", "Gradient"), each = length(second_derivatives)))

# Save Alpha distribution plot
pdf(file = "plots/alpha_fisher_plot.pdf")
ggplot(df_alpha, aes(x = value, color = type, fill = type)) +
  geom_density(alpha = 0.4) +
  #geom_vline(xintercept = alpha_fisher_entry, color = "black", lwd = 1) +
  labs(x = TeX("${I}(\\theta)_{1,1}$"),
       y = "Density") +
  theme_minimal() +
  theme(legend.position = "none")
dev.off()

# MU
mu_fisher_entry <- mu_fisher(3, 2, 0.01)
i <- 2
second_derivatives <- sapply(results$fisher_info_second_derivatives_list, function(m) m[i, i])
gradients <- sapply(results$fisher_info_gradient_list, function(m) m[i, i])

df_mu <- data.frame(value = c(second_derivatives, gradients),
                    type = rep(c("Second Derivative", "Gradient"), each = length(second_derivatives)))

# Save Mu distribution plot
pdf(file = "plots/mu_fisher_plot.pdf")
ggplot(df_mu, aes(x = value, fill = type, color = type)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = mu_fisher_entry, color = "black", lwd = 1) +
  labs(x = TeX("${I}(\\theta)_{2,2}$"),
       y = "Density") +
  theme_minimal() +
  theme(legend.position = "none")
dev.off()

# SIGMA^2
sigma_2_fisher_entry <- sigma2_fisher(3, 2, 0.01)
i <- 3
second_derivatives <- sapply(results$fisher_info_second_derivatives_list, function(m) m[i, i])
gradients <- sapply(results$fisher_info_gradient_list, function(m) m[i, i])

df_sigma2 <- data.frame(value = c(second_derivatives, gradients),
                        type = rep(c("Second Derivative", "Gradient"), each = length(second_derivatives)))

# Save Sigma^2 distribution plot
pdf(file = "plots/sigma2_fisher_plot.pdf")
ggplot(df_sigma2, aes(x = value, color = type, fill = type)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = sigma_2_fisher_entry, color = "black", lwd = 1) +
  labs(x = TeX("${I}(\\theta)_{3,3}$"),
       y = "Density") +
  theme_minimal() +
  theme(legend.position = "none")
dev.off()