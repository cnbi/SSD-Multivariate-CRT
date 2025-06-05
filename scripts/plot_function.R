################ PLOTS ###########################

#my_colour: 
density_plot <- function(data, BF_thresh, save = FALSE, x_axis, limit, my_colour) {
    
    # Libraries
    library(ggplot2)
    library(scales) # For scaling log10
    library(cowplot) #Side-by-side plot+
    
    # Plotting
    density_plot <- ggplot(data = data.frame(data), aes(x = {{x_axis}})) +
        geom_density(adjust = 2, color = "black", size = 0.3) + 
        geom_vline(aes(xintercept = BF_thresh), color = "black", 
                   linetype = "dashed") +
        xlim(0, limit)
    # Density for shade
    shade <- ggplot_build(density_plot)
    x1 <- min(which(shade$data[[1]]$x > (BF_thresh + 0.5)))
    x2 <- max(which(shade$data[[1]]$x <= limit))
    density_plot <- density_plot +
        geom_area(data = data.frame(x = shade$data[[1]]$x[x1:x2], y = shade$data[[1]]$y[x1:x2]),
                  linetype = 1, color = "black",
                  aes(x = x, y = y), fill = my_colour, alpha = 0.6) +
        ylab("Density") +
        theme_classic() +
        xlab("Bayes Factor") +
        # Annotations
        annotate("text", x = 9, y = 0.03, 
                 label = paste(expression(eta), " = ", data$Proportion.BF1c), 
                 color = "black", parse = TRUE, size = 2.8) +
        annotate(geom = "text", x = 3.5, y = 0.025,
                 label = paste(expression(BF[t][h][r][e][s]), " = ", BF_thresh), 
                 angle = 90, color = "black", size = 3.5)
    # Save plot
    if (save == TRUE) {
        ggsave("both_densities.png", width = 6, height = 4, dpi = 600)
        ggsave("both_densities.tiff", width = 6, height = 4, dpi = 600)
    }
    # Return the plot
    return(density_plot)
}


## Saving side-by-side plots with common y axis
plot_grid(H0_true, H1_true + theme(axis.title.y = element_blank()))
ggsave("both_densities.png", width = 6, height = 4, dpi = 600)
ggsave("both_densities.tiff", width = 6, height = 4, dpi = 600)