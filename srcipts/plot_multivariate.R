################### Plot multivariate #############################


# Generate multivariate normal data with 500 samples
set.seed(15042025)  # For reproducibility
n_samples <- 500
data <- mvrnorm(n_samples, mu = estimates, Sigma = sigma)

# Convert the data into a data frame for ggplot2
data_df <- as.data.frame(data)
colnames(data_df) <- c("Y1", "Y2")

# Plot the data
ggplot(data_df, aes(x = Y1, y = Y2)) +
    geom_point(alpha = 0.6) +  # Scatter plot of the points
    stat_ellipse(level = 0.95, aes(color = "Confidence Ellipse")) +  # 95% confidence ellipse
    theme_minimal() +
    labs(title = "Multivariate Normal Distribution",
         x = "Outcome 1 (Y1)", y = "Outcome 2 (Y2)") +
    theme(legend.position = "none")

#3D plot
library(MASS)
library(plotly)
set.seed(123)  # For reproducibility
n_samples <- 500
data <- mvrnorm(n_samples, mu = c(0, 0), Sigma = sigma)

# Convert the data into a data frame
data_df <- as.data.frame(data)
colnames(data_df) <- c("Y1", "Y2")

# Calculate the joint density using dmvnorm
density <- dmvnorm(data_df, mean = estimates, sigma = sigma)

# Add the density as the third dimension
data_df$Z <- density

# Create the 3D scatter plot
plot_ly(data_df, x = ~Y1, y = ~Y2, z = ~Z, 
        type = "scatter3d", mode = "markers", 
        marker = list(size = 3, color = ~Z, colorscale = "Viridis", showscale = TRUE)) %>%
    layout(title = "3D Plot of Two Outcomes with Density as Third Dimension",
           scene = list(xaxis = list(title = "Outcome 1 (Y1)"),
                        yaxis = list(title = "Outcome 2 (Y2)"),
                        zaxis = list(title = "Density")))