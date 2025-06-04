################### Plot multivariate #############################

# Intersection-union test ----------------------------------------
## Simple 2D plot ------------------------------------------------- 
# Libraries
library(MASS) #For multivariate data generation
library(ggplot2)

# Generate multivariate normal data with 500 samples
set.seed(15042025)  # For reproducibility
n_samples <- 500
data <- mvrnorm(n_samples, mu = estimates, Sigma = sigma)
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

## Interactive 3D plot --------------------------------------------
# Libraries
library(MASS)    #For multivariate data generation
library(plotly)
library(mvtnorm) # Multivariate normal density
library(ggExtra)

# Data generation
set.seed(123)  # For reproducibility
n_samples <- 500
data <- mvrnorm(n_samples, mu = c(0, 0), Sigma = sigma)
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

# 2D plot with marginal distributions
# Load required libraries
library(ggplot2)
library(MASS)  # for mvrnorm

# Define parameters for the bivariate normal distribution
mu <- c(0, 0)  # Mean vector for the two outcomes
sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)  # Covariance matrix

# Generate a grid of points
set.seed(2905)
data <- mvrnorm(1000, mu = mu, Sigma = sigma)
df <- data.frame(Outcome1 = data[, 1], Outcome2 = data[, 2])

base_plot <- ggplot(df, aes(x = Outcome1, y = Outcome2)) +
    geom_point(alpha = 0.5, color = "#40B0A6") +
    labs(title = "Joint Distribution with Marginals", 
         x = "Outcome 1 (β1)", y = "Outcome 2 (β2)") +
    coord_cartesian(xlim = c(-3, 3), ylim = c(-3, 3))  # Highlight region where β1, β2 > 0

# Add the shaded region where β1 > 0 and β2 > 0
base_plot <- base_plot +
    geom_point(data = subset(df, Outcome1 > 0 & Outcome2 > 0), aes(x = Outcome1, y = Outcome2),
               color = "#F4B942", alpha = 0.5, size = 2)

final_plot <-  ggMarginal(base_plot, type = "density")
final_plot <-  ggMarginal(base_plot, type = "histogram")

final_plot <- ggMarginal(base_plot, 
           type = "histogram", 
           fill = "#85C0F9",  # Custom fill colour for the marginal histograms
           color = "#6B9AC4",  # Custom border colour for histograms
           size = 5) 

final_plot


# Add marginal histograms using ggExtra package (no additional arguments)
# df$Outcome1_color <- ifelse(df$Outcome1 > 0, "#DEAB44", "#85C0F9")
# df$Outcome2_color <- ifelse(df$Outcome2 > 0, "#DEAB44", "#85C0F9")
# 
# final_plot <- ggMarginal(base_plot, type = "histogram",
#                          fill = df$Outcome1_color,
#                          color = "#6B9AC4",
#                          margins = "x") +
#     ggMarginal(base_plot, type = "histogram",
#                fill = df$Outcome2_color,
#                color = "#6B9AC4",
#                margins = "y")
# 
# final_plot <- base_plot + 
#     ggside::geom_xsidedensity(aes(fill = df$Outcome1_color), 
#                               size = 1, alpha = 0.6, color = "black") + 
#     ggside::geom_ysidedensity(aes(fill = df$Outcome2_color), 
#                               size = 1, alpha = 0.6, color = "black") + 
#     scale_fill_identity() +  # Use the colours directly without scaling
#     theme_minimal() +
#     theme(legend.position = "none")
# 
# # Display the final plot
# final_plot
# 
# 
# # Create a new data frame for each part of the outcome distributions
# final_plot <- base_plot + 
#     ggside::geom_xsidedensity(aes(x = df$Outcome1, fill = ifelse(Outcome1 > 0, "darkblue", "lightblue")),
#                               size = 1, alpha = 0.6, color = "black") + 
#     
#     # Density plot for Outcome 2 (Y-axis), with colour change at the threshold
#     ggside::geom_ysidedensity(aes(x = df$Outcome2, fill = ifelse(Outcome2 > 0, "darkblue", "lightblue")),
#                               size = 1, alpha = 0.6, color = "black") + 
#     
#     scale_fill_identity() +  # Use the colours directly without scaling
#     theme_minimal() +
#     theme(legend.position = "none")
# 
# # Display the final plot
# final_plot

# Homogeneity of effect size test
# Load required libraries
library(ggplot2)
library(MASS)  # for mvrnorm

# Define parameters for the bivariate normal distribution
mu <- c(0, 0)  # Mean vector for the two outcomes
sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)  # Covariance matrix

# Generate a grid of points
set.seed(2905)
data <- mvrnorm(1000, mu = mu, Sigma = sigma)
df <- data.frame(Outcome1 = data[, 1], Outcome2 = data[, 2])

base_plot <- ggplot(df, aes(x = Outcome1, y = Outcome2)) +
    geom_point(alpha = 0.5, color = "#40B0A6") +
    labs(title = "Joint Distribution with Marginals", 
         x = "Outcome 1 (β1)", y = "Outcome 2 (β2)") +
    coord_cartesian(xlim = c(-3, 3), ylim = c(-3, 3))  # Highlight region where β1, β2 > 0

# Add the shaded region where β1 > 0 and β2 > 0
base_plot <- base_plot +
    geom_point(data = subset(df, Outcome1 > 0 & Outcome2 > 0), aes(x = Outcome1, y = Outcome2),
               color = "#F4B942", alpha = 0.5, size = 2)

final_plot <-  ggMarginal(base_plot, type = "density")
final_plot <-  ggMarginal(base_plot, type = "histogram")

final_plot <- ggMarginal(base_plot, 
                         type = "histogram", 
                         fill = "#85C0F9",  # Custom fill colour for the marginal histograms
                         color = "#6B9AC4",  # Custom border colour for histograms
                         size = 5) 

final_plot