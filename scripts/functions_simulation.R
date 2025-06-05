########### SIMULATION FUNCTIONS ###############


# simulation function for second project
# name_results: characters. Name which will be use fo saving results
# name_times: characters. Name which will be use fo saving time results
run_simulation <- function(Row, name_results, name_times, design_matrix, results_folder, seed){
    # Start time
    start_time <- Sys.time()

    # Actual simulation
    ssd_results <- SSD_mult_CRT(test = design_matrix[Row, "test"], 
                                  effect_sizes = c(design_matrix[Row, "eff_size1"], design_matrix[Row, "eff_size2"]), 
                                  n1 = design_matrix[Row, "n1"],
                                  n2 = design_matrix[Row, "n2"], 
                                  ndatasets = 100,
                                  out_specific_ICC = design_matrix[Row, "out_specific_ICC"], 
                                  intersubj_between_outICC = design_matrix[Row, "intersubj_between_outICC"], 
                                  intrasubj_between_outICC = design_matrix[Row, "intrasubj_between_outICC"],
                                  BF_thresh = design_matrix[Row, "BF_thresh"], eta = 0.8, 
                                  fixed = as.character(design_matrix[Row, "fixed"]), max = 1000, 
                                  batch_size = 1000,
                                  Bayes_pack = as.character(design_matrix[Row, "Bayes_pack"]),
                                seed = seed)
    
    # Save results
    end_time <- Sys.time()
    file_name <- file.path(results_folder, paste0(name_results, Row, ".RDS"))
    saveRDS(ssd_results, file = file_name)
    
    # Save running time
    running_time <- as.numeric(difftime(end_time,start_time, units = "mins"))
    time_name <- file.path(results_folder, paste0(name_times, Row, ".RDS"))
    saveRDS(running_time, file = time_name)
    
    # Clean
    rm(ssd_results)
    gc()
}

# Simulation function



# Simulation parallelised
#required_fx: Vector with the objects and functions required to run the simulation

simulation_parallelised <- function(fx, design_matrix, folder, nclusters, parall,
                                    required_fx){
    
    # Packages
    if (!require(parallel)) {install.packages("parallel")}
    if (!require(foreach)) {install.packages("foreach")}
    if (!require(doParallel)) {install.packages("doParallel")}
    if (!require(dplyr)) {install.packages("dplyr")}
    if(!require(future.apply)){install.packages("future.apply")}
    nrow_design <- nrow(design_matrix)
    
    if (parall == "doParallel"){
        # Detect
        ncluster <- detectCores() / 2
        # Create clusters and register them
        cl <- makeCluster(ncluster)
        registerDoParallel(cl)
        # Export libraries, functions and variables
        clusterExport(cl, required_fx)
        # Distribute rows into clusters
        dentify_clusters <- cut(seq(nrow_design), breaks = ncluster, labels = FALSE)
        rows_divided <- split(seq(nrow_design), 1:ncluster)
        rows_divided <- lapply(seq(nrow_design), function(x) list(x))
        # Apply simulation
        clusterApply(cl, rows_divided, run_simulation)
        # Stop parallelisation
        stopCluster(cl)
        stopImplicitCluster()
        
        
    } else if (parall == "forEach") {
        # Detect
        ncluster <- detectCores() / 2
        # Create clusters and register them
        cl <- makeCluster(ncluster)
        registerDoParallel(cl)
        # Distribute rows into clusters
        rows_divided <- split(seq(nrow_design), 1:ncluster)
        # Export libraries, functions and variables
        clusterExport(cl, required_fx)
        
        # Parallelisation
        foreach(Rows = rows_divided) %dopar% {
            run_simulation(Rows)
        }
        # Stop clusters
        stopCluster(cl)
        stopImplicitCluster()
        
    } else if (parall == "future") {
        plan(multisession)
        future_lapply()
        
    }
    
}

# Collect results
# Collect results in a matrix---
collect_results <- function(design_matrix, results_folder, finding, pair, name_results) {
    rows <-  seq(nrow(design_matrix))
    if (finding == "N2") {
        results_name <- ifelse(missing(name_results), "/ResultsN2Row", paste0("/", name_results))
        file_name <- "final_results_findN2"
    } else if (finding == "N1") {
        results_name <- "/ResultsN1Row"
        file_name <- "final_results_findN1"
    }
    
    new_matrix <- matrix(NA, ncol = 6, nrow = nrow(design_matrix))
    for (row_design in rows) {
        stored_result <- readRDS(paste0(results_folder, results_name, row_design, ".RDS"))
        median.BF1c <- median(stored_result[[6]][, "BF.1c"])
        median.BF1u <- median(stored_result[[6]][, "BF.1u"])
        mean.PMP1c <- mean(stored_result[[6]][, "PMP.1c"])
        eta.BF1c <- stored_result$Proportion.BF1c
        n2 <- stored_result$n2
        n1 <- stored_result$n1
        new_matrix[row_design, ] <- c(median.BF1c, median.BF1u,
                                      mean.PMP1c, eta.BF1c,
                                      n2, n1)
    }
    
    new_matrix <- as.data.frame(cbind(design_matrix, new_matrix))
    colnames(new_matrix) <- c(names(design_matrix), "median.BF1c",
                              "median.BF1c", "mean.PMP1c",
                              "eta.BF1c", "n2.final", "n1.final")
    saveRDS(new_matrix, file = file.path(results_folder, paste0(file_name, ".RDS")))
    
    return(new_matrix)
}

# Collect times in a matrix ----
collect_times <- function(design_matrix, rows = 0, pair, finding, results_folder) {
    rows <- ifelse(rows = 0, seq(nrow(design_matrix)), rows)
    new_matrix <- matrix(NA, nrow = nrow(design_matrix), ncol = 1)
    if (finding == "N2") {
        results_name <- "/timeN2Row"
        file_name <- "final_times_findN2"
    } else if (finding == "N1") {
        results_name <- "/timeN1Row"
        file_name <- "final_times_findN1"
    }
    
    for (row_result in rows) {
        stored_result <- readRDS(paste0(results_folder, results_name, row_result, ".RDS"))
        new_matrix[row_result, ] <- stored_result
    }
    new_matrix <- as.data.frame(cbind(design_matrix, new_matrix))
    colnames(new_matrix) <- c(names(design_matrix), "total.time")
    saveRDS(new_matrix, file = file.path(results_folder, paste0(file_name, type)))
    
    return(new_matrix)
}