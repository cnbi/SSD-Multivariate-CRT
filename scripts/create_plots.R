################## Function for plot ###############################

# Density plot



density_plot <- function(data){
    if (!require("ggplot2")) {install.packages("ggplot2")} #Plots
    if (!require("cowplot")) {install.packages("cowplot")} #Side-by-side plot
    
    data_plot <- data$data_H1
    BF_threshold <- data$BF_thres
    emp_power <- data$eta
    
    H1_true <- ggplot(data = data.frame(data_plot), aes(x = BF.1c)) +
        geom_density(adjust = 2, color = "black", size = 0.3) +
        geom_vline(aes(xintercept = BF_threshold),
                   color = "black", linetype = "dashed") +
        xlim(0, 40)
    shade <- ggplot_build(H1_true)
    x1 <- min(which(shade$data[[1]]$x > BF_threshold+ .05))
    x2 <- max(which(shade$data[[1]]$x <= 40))
    H1_true <- H1_true + 
        geom_area(data = data.frame(x = shade$data[[1]]$x[x1:x2], y = shade$data[[1]]$y[x1:x2]),
                  linetype = 1, color = "black",
                  aes(x = x, y = y), fill = "#D3D3D3", alpha = 0.6) +
        ylab("Density") +
        theme_classic() +
        annotate("text", x = 10, y = 0.02, label = expression(eta == 0.82), color = "black", parse = TRUE, size = 2.8) +
        xlab("Bayes Factor") +
        annotate(geom = "text", x = 3.5, y = 0.0165,label = expression(BF[t][h][r][e][s] == 5), angle = 90, color = "black", size = 3.5)
    

    
    # Other option
    H1_true <- ggplot(data = data.frame(df_H1), aes(x = BF.10)) +
        geom_density(adjust = 2, color = "black", size = 0.3) +  # Density plot
        geom_vline(aes(xintercept = BF_threshold),  # Dashed vertical line
                   color = "black", linetype = "dashed") +
        geom_area(data = subset(data.frame(df_H1), BF.10 >= BF_threshold),  # Shade area to the right of the line
                  aes(x = BF.10, y = ..density..), fill = "#D3D3D3", alpha = 0.6) +
        xlim(0, 40) +
        ylab("Density") +
        xlab("Bayes Factor") +
        theme_classic() +
        annotate("text", x = 10, y = 0.02, label = expression(eta == 0.82), color = "black", parse = TRUE, size = 2.8) +
        annotate(geom = "text", x = 3.5, y = 0.0165, label = expression(BF[t][h][r][e][s] == 5), angle = 90, color = "black", size = 3.5)
    
    # Display plot
    print(H1_true)
    
}