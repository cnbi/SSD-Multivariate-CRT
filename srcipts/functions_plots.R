############# PLOTS ##############################

library(ggplot2)
library(rlang)

# Plot for main effects of factors -----------------------
main_eff_plot <- function(data, x, y, grouping, fixed, val_fixed,
                          color_labs, xlab, ylab, legend_pos) {
    conditions <- paste(fixed, "==", val_fixed, collapse = " & ")
    data_plots <- subset(data, eval(parse(text = conditions)))
    x <- sym(x)
    y <- sym(y)
    grouping <- sym(grouping)
    max_y <- max(data_plots[[as.character(y)]])
    ggplot(data_plots, aes(x = !!x, y = !!y,
                           color = as.factor(!!grouping), 
                           shape = as.factor(!!grouping))) +
        geom_point() + geom_line() + scale_color_brewer(palette = "Set2") +
        scale_fill_brewer("Set2") + labs(color = color_labs, shape = color_labs) +
        xlab(xlab) + ylab(ylab) + theme(legend.position = legend_pos) +
        ylim(0, (max_y + 5))
}

# Plot of interactions
inter_eff_plot <- function(data, x, y, grouping, fixed, val_fixed,
                           color_labs, xlab, ylab, legend_pos, panels,
                           labels_grid) {
    conditions <- paste(fixed, "==", val_fixed, collapse = " & ")
    data_plots <- subset(data, eval(parse(text = conditions)))
    x <- sym(x)
    y <- sym(y)
    grouping <- sym(grouping)
    max_y <- max(data_plots[[as.character(y)]])
    base_plot <- ggplot(data_plots, aes(x = !!x, y = !!y,
                                        color = as.factor(!!grouping), 
                                        shape = as.factor(!!grouping))) +
        geom_point() + geom_line() + scale_color_brewer(palette = "Set2") +
        scale_fill_brewer("Set2") + labs(color = color_labs, shape = color_labs) +
        xlab(xlab) + ylab(ylab) + theme(legend.position = legend_pos) + ylim(0, (max_y + 5))
    if (length(panels) > 1) {
        panel_rows <- sym(panels[1])
        panel_cols <- sym(panels[2])
        if (missing(labels_grid)) {
            final_plot <- base_plot + facet_grid(rows = vars(!!panel_rows), 
                                                 cols = vars(!!panel_cols))
        } else {
            final_plot <- base_plot + facet_grid(rows = vars(!!panel_rows), 
                                                 cols = vars(!!panel_cols),
                                                 labeller = labels_grid)
        }
    } else if (length(panels) == 1) {
        panel_rows <- sym(panels)
        if (missing(labels_grid)) {
            final_plot <- base_plot + facet_grid(rows = vars(!!panel_rows))
        } else {
            final_plot <- base_plot + facet_grid(rows = vars(!!panel_rows),
                                                 labeller = labels_grid)
            
        } 
        return(final_plot)
    }
}

inter_eff_plot(data = final_results_findN2,
               x = "n1.final", y = "n2.final",
               grouping = "BF_thresh",
               fixed = c("eff_size1", "eff_size2",
                         "intersubj_between_outICC",
                         "intrasubj_between_outICC"),
               val_fixed = c(0.3, 0.5, 0.005, 0.2),
               color_labs = "Bayes Factor \nThreshold",
               xlab = "Cluster size", ylab = "Number of clusters",
               legend_pos = "right", panels = c("out_specific_ICC"))


# TODO: Add plotly or ggiraph option