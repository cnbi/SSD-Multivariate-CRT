###################### SIMULATION ########################################333

# Simulation for intersection-union test
# Design matrix
# General factors
eff_size <- c(0.3, 0.5, 0.7, 0.9)
outcome_icc <- c(0.01, 0.025, 0.05)
intersubj_icc <- c(0.005, 0.025)
intrasubj_icc <- c(0.2, 0.5)
bf_thresh <- c(5, 10)
eta <- c(0.8)
# test
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
design_matrix_n2 <- expand.grid(outcome_icc, intersubj_icc, intrasubj_icc, bf_thresh,
                                 eta, fixed, n1, n2, test, bf_pack)
effect_sizes  <- matrix(rep(t(effect_sizes), nrow(design_matrix_n2)),ncol = ncol(effect_sizes), byrow = TRUE)
effect_sizes <- effect_sizes[order(effect_sizes[, 1], effect_sizes[, 2]), ]
design_matrix_n2 <- cbind.data.frame(effect_sizes, design_matrix_n2)
colnames(design_matrix_n2) <- c("eff_size1", "eff_size2", "out_specific_ICC", "intersubj_between_outICC", "intrasubj_between_outICC", "BF_thresh",
                                "eta", "fixed", "n1", "n2", "test", "Bayes_pack")

# Finding cluster size
n2 <- c(6, 10, 40)
fixed <- c("n2")
n1 <- 15
design_matrix_n1 <- expand.grid(eff_size, outcome_icc, intersubj_icc, intrasubj_icc, bf_thresh,
                                 eta, fixed, n2, n1, test, bf_pack)
colnames(design_matrix_n1) <- c("effect_sizes", "out_specific_ICC", "intersubj_between_outICC", "intrasubj_between_outICC", "BF_thresh",
                                "eta", "fixed", "n2", "n1", "test", "Bayes_pack")

# Source
source("scripts/functions_simulation.R")
source("scripts/find_n2_multiv.R")

#Libraries
library(dplyr)   # Format tables
library(purrr)    # Format tables

# Creation folder for results
folder_results <- "srcipts/findN2_iu"
if (!dir.exists(folder_results)) {dir.create(folder_results)}

# Loop for every row
for (Row in 1) {
    run_simulation(Row, name_results = arg_fx[1], name_times = arg_fx[2], 
                          design_matrix = design_matrix_n2, results_folder = folder, seed = 2106)
}

# Collect results
collect_results(design_matrix = design_matrix_n2, results_folder = folder_results, finding = "N2", 
                name_results = "findN2_iu_")
index_missingBF <- which(is.na(final_results_findN2$median.BF1c))
missing_BF <- final_results_findN2[which(is.na(final_results_findN2$median.BF1c)), ]
# Loop for every row
for (Row in index_missingBF) {
    seed <- 2106
    run_simulation(Row, name_results = "findN2_iu_", name_times = "time_findN2_iu", 
                   design_matrix = design_matrix_n2, results_folder = folder_results, seed)
}

# Parallelise version with future
# Source
source("scripts/functions_simulation.R")
source("scripts/find_n2_multiv.R")
# Creation folder for results
folder_results <- "scripts/FindN2_iu_new"
if (!dir.exists(folder_results)) {dir.create(folder_results)}
# Run simulation
arg_fx <- c("FindN2_IU_", "TimeN2_IU", 2106)
simulation_parallelised(design_matrix = design_matrix_n2, folder = folder_results, nclusters = 5,
                        parall = "do", required_fx = arg_fx)

simulation_parallelised(design_matrix = design_matrix_n2, folder = folder_results, nclusters = 5,
                        parall = "forEach", required_fx = arg_fx)


# Detect -----------------------------------------------
nclusters <- 5
if (missing(nclusters)) {
    ncluster <- detectCores() / 2
} else {
    ncluster <- nclusters
}
nrow_design <- design_matrix_n2
required_fx <- arg_fx
design_matrix <- design_matrix_n2
folder <- folder_results
# Create clusters and register them
cl <- makeCluster(ncluster)
registerDoParallel(cl)
# Distribute rows into clusters
rows_divided <- split(seq(nrow_design), 1:ncluster)
# Export libraries, functions and variables
# Parallelisation
future::plan(multisession, workers = 4)
rows_to_run <- 1:nrow(design_matrix)
library(doFuture)
foreach(rows_to_run, .export = c("run_simulation", "design_matrix", "required_fx", "folder", "SSD_mult_CRT")) %dofuture% {
    run_simulation(Rows, name_results = required_fx[1], name_times = required_fx[2],
                   design_matrix = design_matrix, results_folder = folder,
                   seed = as.integer(required_fx[3]))
}
# Stop clusters
stopCluster(cl)
stopImplicitCluster()

# Future------------------
if (missing(nclusters)) {
    nclusters <- detectCores() / 2
} else {
    nclusters <- nclusters
}
plan(multisession, workers = nclusters)
rows_to_run <- 1:nrow(design_matrix)
future_lapply(rows_to_run, function(Row)
    run_simulation(Row = Row, 
                   name_results = required_fx[1], 
                   name_times = required_fx[2],
                   design_matrix = design_matrix,
                   results_folder = folder,
                   seed = as.integer(required_fx[3])
    )
)
    
future_lapply(rows_to_run, run_simulation, 
              Row = rows_to_run, 
              name_results = required_fx[1], 
              name_times = required_fx[2],
              design_matrix = design_matrix,
              results_folder = folder,
              seed = as.integer(required_fx[3]))
