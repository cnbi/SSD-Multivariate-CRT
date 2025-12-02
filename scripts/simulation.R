###################### SIMULATION ########################################333

# Simulation for intersection-union test------------------------------------------
## Design matrix-----------------------------------------------------------------
# General factors
eff_size <- c(0.3, 0.5, 0.7, 0.9)
outcome_icc <- c(0.01, 0.05)
intersubj_icc <- c(0.005, 0.025)
intrasubj_icc <- c(0.2, 0.5)
pmp_thresh <- c(0.9, 0.95)
eta <- c(0.8)

# Test-specific factors
test <- c("intersection-union")
effect_sizes <- matrix(NA, nrow = 6, ncol = 2)
i <- 1
if (test == "intersection-union") {
    for (a in 1:length(eff_size)) {
        eff1 <- eff_size[a]
        for (b in 2:length(eff_size)) {
            eff2 <- eff_size[b]
            if (!eff1 == eff2) {
                if (eff1 < eff2) {
                    effect_sizes[i, ] <- c(eff1, eff2)
                    i <- i + 1
                }
            } 
        }
        
    }
}
bf_pack <- c("bain")

# Finding number of clusters
n1 <- c(5, 15, 30) 
n2 <- 30
fixed <- c("n1")
design_matrix_n2 <- expand.grid(outcome_icc, intersubj_icc, intrasubj_icc, pmp_thresh,
                                eta, fixed, n1, n2, test, bf_pack)
effect_sizes  <- matrix(rep(t(effect_sizes), nrow(design_matrix_n2)),ncol = ncol(effect_sizes), byrow = TRUE)
effect_sizes <- effect_sizes[order(effect_sizes[, 1], effect_sizes[, 2]), ]
design_matrix_n2 <- cbind.data.frame(effect_sizes, design_matrix_n2)
colnames(design_matrix_n2) <- c("eff_size1", "eff_size2", "out_specific_ICC", "intersubj_between_outICC", "intrasubj_between_outICC", "pmp_thresh",
                                "eta", "fixed", "n1", "n2", "test", "Bayes_pack")

# Finding cluster size
n2 <- c(10, 40, 60) # Total number of clusters
fixed <- c("n2")
n1 <- 15
design_matrix_n1 <- expand.grid(eff_size, outcome_icc, intersubj_icc, intrasubj_icc, pmp_thresh,
                                eta, fixed, n2, n1, test, bf_pack)
colnames(design_matrix_n1) <- c("effect_sizes", "out_specific_ICC", "intersubj_between_outICC", "intrasubj_between_outICC", "pmp_thresh",
                                "eta", "fixed", "n2", "n1", "test", "Bayes_pack")

# Source
source("functions_simulation.R")
source("find_n2_multiv.R")

#Libraries
library(dplyr)   # Format tables
library(purrr)    # Format tables

# Creation folder for results
folder_results <- "experiments"
if (!dir.exists(folder_results)) {dir.create(folder_results)}
arg_fx <- c("FindN2_IU_", "TimeN2_IU", 610)

## Loop for every row-----------------------------------------------------------
for (Row in 3) {
    run_simulation(Row, name_results = arg_fx[1], name_times = arg_fx[2], 
                   design_matrix = design_matrix_n2, results_folder = folder_results, seed = 17)
}

## Collect results------------------------------------------------------------
collect_results(design_matrix = design_matrix_n2, results_folder = folder_results, finding = "N2", 
                name_results = "findN2_iu_")
index_missingBF <- which(is.na(final_results_findN2$median.BF1c))
missing_BF <- final_results_findN2[which(is.na(final_results_findN2$median.BF1c)), ]
# Loop for every row
for (Row in 81) {
    seed <- 2106
    options(error = recover)
    run_simulation(Row, name_results = "FindN2_IU_", name_times = "TimeN2_IU", 
                   design_matrix = design_matrix_n2, results_folder = folder_results, seed)
}

## Parallelise version with future-----------------------------------------------
# Source
source("functions_simulation.R")
source("find_n2_multiv.R")

# Creation folder for results
folder_results <- "FindN2_iu_new"
if (!dir.exists(folder_results)) {dir.create(folder_results)}
# Run simulation

arg_fx <- c("FindN2_IU_", "TimeN2_IU_", 810) #Name of results, name of time, seed

simulation_parallelised(design_matrix = design_matrix_n2[1:5, ], folder = folder_results, nclusters = 5,
                        parall = "future", required_fx = arg_fx)


## Collect results in one matrix --------------------------------------------------
results_iu <- collect_results(design_matrix_n2, finding = "N2", results_folder = "IU", 
                name_results = "FindN2_IU_", test = test, file_name = "results_FindN2_IU")

times_iu <- collect_times(design_matrix = design_matrix_n2, results_folder = "IU", 
                          finding = "N2", name_results = "TimeN2_IU", test = test,
                          file_name = "times_FindN2_IU")

 # Simulation fo Homogeneity of effect sizes-------------------------------------
# General factors
outcome_icc <- c(0.01, 0.05)
intersubj_icc <- c(0.005, 0.025)
intrasubj_icc <- c(0.2, 0.5)
pmp_thresh <- c(0.9, 0.95)
eta <- c(0.8)

# Test-specific factors
test <- c("homogeneity")
delta <- c(0.2, 0.3)
eff_size <- c(0.3, 0.6, 0.9)
effect_sizes <- matrix(NA, nrow = 3, ncol = 2)
i <- 1
for (a in 1:length(eff_size)) {
    eff1 <- eff_size[a]
    eff2 <- eff1 - 0.1
    effect_sizes[i, ] <- c(eff1, eff2)
    i <- i + 1
}

bf_pack <- c("bain")
# Finding number of clusters
n1 <- c(5, 15, 30) 
n2 <- 30
fixed <- c("n1")
design_matrix_n2 <- expand.grid(delta, outcome_icc, intersubj_icc, intrasubj_icc, pmp_thresh,
                                eta, fixed, n1, n2, test, bf_pack)
effect_sizes  <- matrix(rep(t(effect_sizes), nrow(design_matrix_n2)),ncol = ncol(effect_sizes), byrow = TRUE)
effect_sizes <- effect_sizes[order(effect_sizes[, 1], effect_sizes[, 2]), ]
design_matrix_n2 <- cbind.data.frame(effect_sizes, design_matrix_n2)
colnames(design_matrix_n2) <- c("eff_size1", "eff_size2", "delta", "out_specific_ICC", "intersubj_between_outICC", "intrasubj_between_outICC", "pmp_thresh",
                                "eta", "fixed", "n1", "n2", "test", "Bayes_pack")

# Finding cluster size
n2 <- c(10, 40, 60) # Total number of clusters
fixed <- c("n2")
n1 <- 15
design_matrix_n1 <- expand.grid(eff_size, outcome_icc, intersubj_icc, intrasubj_icc, pmp_thresh,
                                eta, fixed, n2, n1, test, bf_pack)
colnames(design_matrix_n1) <- c("effect_sizes", "out_specific_ICC", "intersubj_between_outICC", "intrasubj_between_outICC", "pmp_thresh",
                                "eta", "fixed", "n2", "n1", "test", "Bayes_pack")

# Source
source("functions_simulation.R")
source("find_n2_multiv.R")

#Libraries
library(dplyr)   # Format tables
library(purrr)    # Format tables

# Creation folder for results
folder_results <- "homogeneity"
if (!dir.exists(folder_results)) {dir.create(folder_results)}
arg_fx <- c("FindN2_homog_", "TimeN2_homog_", 610)

# Change name
change <- list.files(folder_results, recursive = TRUE, full.names = TRUE)
change <- change[252:502]
change <- str_sort(change , numeric = TRUE)
index <- c(252:502)
change_df <- cbind(index, change)
change_df <- change_df[order(as.numeric(change_df[ , 1]), decreasing = TRUE), ]
change_vector <- change_df[, 2]
for (file in change_vector) {
    file_name <- basename(file)
    file_name_alone <- gsub(".RDS$", "", file_name)
    number <- gsub("TimeN2_homoge_", "", file_name_alone)
    actual_number <- as.integer(number) + 36
    new_name <- paste0(folder_results, "/TimeN2_homoge_", actual_number, ".RDS")
    file.rename(file, new_name)
}
## Parallelise version with future-----------------------------------------------
# Source
source("functions_simulation.R")
source("find_n2_multiv.R")

# Creation folder for results
folder_results <- "FindN2_iu_new"
if (!dir.exists(folder_results)) {dir.create(folder_results)}
# Run simulation

arg_fx <- c("FindN2_IU_", "TimeN2_IU_", 810) #Name of results, name of time, seed
simulation_parallelised(design_matrix = design_matrix_n2, folder = folder_results, nclusters = 8,
                        parall = "Parallel", required_fx = arg_fx)

simulation_parallelised(design_matrix = design_matrix_n2[1:5, ], folder = folder_results, nclusters = 5,
                        parall = "forEach", required_fx = arg_fx)

simulation_parallelised(design_matrix = design_matrix_n2[1:5, ], folder = folder_results, nclusters = 5,
                        parall = "future", required_fx = arg_fx)


# Simulation for omnibus test-----------------------------------------------
# General factors
outcome_icc <- c(0.01, 0.05)
intersubj_icc <- c(0.005, 0.025)
intrasubj_icc <- c(0.2, 0.5)
pmp_thresh <- c(0.9, 0.95)
eta <- c(0.8)

# Test-specific factors
test <- c("omnibus")
eff_size <- c(0.2, 0.3, 0.5, 0.7, 0.9)
effect_sizes <- matrix(NA, nrow = 10, ncol = 2)
i <- 1
for (a in 1:length(eff_size)) {
    eff1 <- eff_size[a]
    for (b in 2:length(eff_size)) {
        eff2 <- eff_size[b]
        if (!eff1 == eff2) {
            if (eff1 < eff2) {
                effect_sizes[i, ] <- c(eff1, eff2)
                i <- i + 1
            }
        } 
    }
}

bf_pack <- c("bain")

# Finding number of clusters
n1 <- c(5, 15, 30) 
n2 <- 30
fixed <- c("n1")
design_matrix_n2 <- expand.grid(outcome_icc, intersubj_icc, intrasubj_icc, pmp_thresh,
                                eta, fixed, n1, n2, test, bf_pack)
effect_sizes  <- matrix(rep(t(effect_sizes), nrow(design_matrix_n2)),ncol = ncol(effect_sizes), byrow = TRUE)
effect_sizes <- effect_sizes[order(effect_sizes[, 1], effect_sizes[, 2]), ]
design_matrix_n2 <- cbind.data.frame(effect_sizes, design_matrix_n2)
colnames(design_matrix_n2) <- c("eff_size1", "eff_size2", "out_specific_ICC", "intersubj_between_outICC", "intrasubj_between_outICC", "pmp_thresh",
                                "eta", "fixed", "n1", "n2", "test", "Bayes_pack")

# Finding cluster size
n2 <- c(10, 40, 60) # Total number of clusters
fixed <- c("n2")
n1 <- 15
design_matrix_n1 <- expand.grid(eff_size, outcome_icc, intersubj_icc, intrasubj_icc, pmp_thresh,
                                eta, fixed, n2, n1, test, bf_pack)
colnames(design_matrix_n1) <- c("effect_sizes", "out_specific_ICC", "intersubj_between_outICC", "intrasubj_between_outICC", "pmp_thresh",
                                "eta", "fixed", "n2", "n1", "test", "Bayes_pack")

# Source
source("functions_simulation.R")
source("find_n2_multiv.R")

#Libraries
library(dplyr)   # Format tables
library(purrr)    # Format tables

# Creation folder for results
folder_results <- "omnibus"
if (!dir.exists(folder_results)) {dir.create(folder_results)}
arg_fx <- c("FindN2_omni_", "TimeN2_omni_", 810)

# Source
source("functions_simulation.R")

## Parallelised simulation
simulation_parallelised(design_matrix = design_matrix_n2, folder = folder_results, nclusters = 160,
                        parall = "future", required_fx = arg_fx) # Name of results file, Name time files, master.seed

## Collect results
results_omni <- collect_results(design_matrix = design_matrix_n2, results_folder = "omnibus", finding = "N2", 
                name_results = "FindN2_omni_", test = "omnibus", file_name = "results_FindN2_omni")
times_omni <- collect_times(design_matrix = design_matrix_n2, results_folder = "omnibus",
                            finding = "N2", name_results = "TimeN2_omni_",
                            test = "omnibus", file_name = "times_FindN2_omni")


# Plots-------------------------------------------------------------------------
library(ggplot2)
#ICCs
results_iu <- results_FindN2_IU
rho_labs_intra <- c("Intrasubject between-out: 0.2", "Intrasubject between-out: 0.5")
names(rho_labs_intra) <- c("0.2", "0.5")
rho_labs_inter <- c("Intersubject between-out: 0.005", "Intersubject between-out: 0.025")
names(rho_labs_inter) <- c("0.005", "0.025")
results_iu_plot <- results_iu[(results_iu$eff_size1 == 0.3) & 
                                  (results_iu$eff_size2 == 0.7) & (results_iu$pmp_thresh==0.95) , ]
base <- ggplot(results_iu_plot, aes(x = n1.final, y = n2.final,
                       color = as.factor(out_specific_ICC), 
                       shape = as.factor(out_specific_ICC))) +
    geom_point() + geom_line() + scale_color_brewer(palette = "Set2") +
    scale_fill_brewer("Set2") + labs(color = "Outcome-specific", shape = "Outcome-specific") +
    xlab("Cluster size") + ylab("Number of clusters") + 
    theme(legend.position = "bottom") + ylim(0, (90 + 5))

plot_iccs <- base + facet_grid(rows = vars(intersubj_between_outICC ), 
                       cols = vars(intrasubj_between_outICC),
                  labeller = labeller(intersubj_between_outICC = rho_labs_inter,
                                      intrasubj_between_outICC = rho_labs_intra))

# Save plot with proportions for slide
ggsave(plot = plot_iccs, filename = "plot_iccs.png", width = 13.33,
       height = 7.5,
       dpi = 300)

# Effect sizes
eff_size1_lab <- c("Effect size 1: 0.3", "Effect size 1: 0.5", "Effect size 1:0.7")
names(eff_size1_lab) <- c("0.3", "0.5", "0.7")
eff_size2_lab <- c("Effect size 2: 0.5", "Effect size 2: 0.7", "Effect size 2: 0.9")
names(eff_size2_lab) <- c("0.5", "0.7", "0.9")
results_iu_plot <- results_iu[(results_iu$intersubj_between_outICC==0.025 ) & 
                                  (results_iu$intrasubj_between_outICC==0.5) & (results_iu$pmp_thresh==0.95) , ]
base <- ggplot(results_iu_plot, aes(x = n1.final, y = n2.final,
                                    color = as.factor(out_specific_ICC), 
                                    shape = as.factor(out_specific_ICC))) +
    geom_point() + geom_line() + scale_color_brewer(palette = "Set2") +
    scale_fill_brewer("Set2") + labs(color = "Outcome-specific", shape = "Outcome-specific") +
    xlab("Cluster size") + ylab("Number of clusters") + 
    theme(legend.position = "bottom") + ylim(0, (90 + 5))

plot_effsiz <- base + facet_grid(rows = vars(eff_size1), 
                  cols = vars(eff_size2),
                  labeller = labeller(eff_size1 = eff_size1_lab,
                                      eff_size2 = eff_size2_lab))

# Save plot with proportions for slide
ggsave(plot = plot_effsiz, filename = "plot_effsiz.png", width = 13.33,
       height = 7.5,
       dpi = 300)

# Thresholds
results_iu_plot <- results_iu[(results_iu$intersubj_between_outICC==0.025 ) & 
                                  (results_iu$intrasubj_between_outICC==0.5) & (results_iu$out_specific_ICC==0.05) , ]
results_iu_plot <- results_iu_plot[results_iu_plot$eff_size2 == 0.9, ]
base <- ggplot(results_iu_plot, aes(x = n1.final, y = n2.final,
                                    color = as.factor(pmp_thresh), 
                                    shape = as.factor(pmp_thresh))) +
    geom_point() + geom_line() + scale_color_brewer(palette = "Set2") +
    scale_fill_brewer("Set2") + labs(color = "PMP threshold", shape = "PMP threshold") +
    xlab("Cluster size") + ylab("Number of clusters") + 
    theme(legend.position = "bottom") + ylim(0, (90 + 5))

plot_thres <- base + facet_grid(cols = vars(eff_size1),
                  labeller = labeller(eff_size1 = eff_size1_lab,
                                      eff_size2 = eff_size2_lab))
# Save plot with proportions for slide
ggsave(plot = plot_thres, filename = "plot_thres.png", width = 13.33,
       height = 7.5,
       dpi = 300)

# OMNIBUS #
#ICCs
results_omni <- results_FindN2_omni
rho_labs_intra <- c("Intrasubject between-out: 0.2", "Intrasubject between-out: 0.5")
names(rho_labs_intra) <- c("0.2", "0.5")
rho_labs_inter <- c("Intersubject between-out: 0.005", "Intersubject between-out: 0.025")
names(rho_labs_inter) <- c("0.005", "0.025")
results_omni_plot <- results_omni[(results_omni$eff_size1==0.3) & 
                                  (results_omni$eff_size2==0.7) & (results_omni$pmp_thresh==0.95) , ]
base <- ggplot(results_omni_plot, aes(x = n1.final, y = n2.final,
                                    color = as.factor(out_specific_ICC), 
                                    shape = as.factor(out_specific_ICC))) +
    geom_point() + geom_line() + scale_color_brewer(palette = "Set2") +
    scale_fill_brewer("Set2") + labs(color = "Outcome-specific", shape = "Outcome-specific") +
    xlab("Cluster size") + ylab("Number of clusters") + 
    theme(legend.position = "bottom") + ylim(0, (90 + 5))

plot_iccs_omni <- base + facet_grid(rows = vars(intersubj_between_outICC ), 
                  cols = vars(intrasubj_between_outICC),
                  labeller = labeller(intersubj_between_outICC = rho_labs_inter,
                                      intrasubj_between_outICC = rho_labs_intra))
# Save plot with proportions for slide
ggsave(plot = plot_iccs_omni, filename = "plot_iccs_omni.png", width = 13.33,
       height = 7.5,
       dpi = 300)

# Effect sizes
eff_size1_lab <- c("Effect size 1: 0.2","Effect size 1: 0.3", "Effect size 1: 0.5",
                   "Effect size 1: 0.7")
names(eff_size1_lab) <- c("0.2", "0.3", "0.5", "0.7")
eff_size2_lab <- c("Effect size 2: 0.3", "Effect size 2: 0.5", "Effect size 2: 0.7", "Effect size 2: 0.9")
names(eff_size2_lab) <- c("0.3", "0.5", "0.7", "0.9")
results_omni_plot <- results_omni[(results_omni$intersubj_between_outICC==0.025 ) & 
                                  (results_omni$intrasubj_between_outICC==0.5) &
                                    (results_omni$pmp_thresh==0.95) , ]
base <- ggplot(results_omni_plot, aes(x = n1.final, y = n2.final,
                                    color = as.factor(out_specific_ICC), 
                                    shape = as.factor(out_specific_ICC))) +
    geom_point() + geom_line() + scale_color_brewer(palette = "Set2") +
    scale_fill_brewer("Set2") + labs(color = "Outcome-specific", shape = "Outcome-specific") +
    xlab("Cluster size") + ylab("Number of clusters") + 
    theme(legend.position = "bottom") + ylim(0, (250 + 5))

plot_effsiz_omni <- base + facet_grid(rows = vars(eff_size1), 
                  cols = vars(eff_size2),
                  labeller = labeller(eff_size1 = eff_size1_lab,
                                      eff_size2 = eff_size2_lab))
# Save plot with proportions for slide
ggsave(plot = plot_effsiz_omni, filename = "plot_effsiz_omni.png", width = 13.33,
       height = 7.5,
       dpi = 300)

# Thresholds
results_omni_plot <- results_omni[(results_omni$intersubj_between_outICC==0.025 ) & 
                                  (results_omni$intrasubj_between_outICC==0.5) & (results_omni$out_specific_ICC==0.05) , ]
results_omni_plot <- results_omni_plot[results_omni_plot$eff_size2 == 0.9, ]
base <- ggplot(results_omni_plot, aes(x = n1.final, y = n2.final,
                                    color = as.factor(pmp_thresh), 
                                    shape = as.factor(pmp_thresh))) +
    geom_point() + geom_line() + scale_color_brewer(palette = "Set2") +
    scale_fill_brewer("Set2") + labs(color = "PMP threshold", shape = "PMP threshold") +
    xlab("Cluster size") + ylab("Number of clusters") + 
    theme(legend.position = "bottom") + ylim(0, (200 + 5))

plot_thres_omni <- base + facet_grid(rows = vars(eff_size2), 
                  cols = vars(eff_size1),
                  labeller = labeller(eff_size1 = eff_size1_lab,
                                      eff_size2 = eff_size2_lab))
# Save plot with proportions for slide
ggsave(plot = plot_thres_omni, filename = "plot_thres_omni.png", width = 13.33,
       height = 7.5,
       dpi = 300)

#HOMOGENEITY
