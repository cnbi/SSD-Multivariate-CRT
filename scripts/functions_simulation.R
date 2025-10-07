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
                                ndatasets = 500,
                                out_specific_ICC = design_matrix[Row, "out_specific_ICC"], 
                                intersubj_between_outICC = design_matrix[Row, "intersubj_between_outICC"], 
                                intrasubj_between_outICC = design_matrix[Row, "intrasubj_between_outICC"],
                                pmp_thresh = design_matrix[Row, "pmp_thresh"], eta = design_matrix[Row, "eta"], 
                                fixed = as.character(design_matrix[Row, "fixed"]), max = 300, 
                                Bayes_pack = as.character(design_matrix[Row, "Bayes_pack"]),
                                master.seed = as.integer(seed))
    
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
#required_fx: Vector with the objects and functions required to run the simulation.
## It includes name_results, name_times, seed, function needed.

simulation_parallelised <- function(design_matrix, folder, nclusters, parall,
                                    required_fx){
    
    # Packages
    if (!require(parallel)) {install.packages("parallel")}
    if (!require(foreach)) {install.packages("foreach")}
    if (!require(doParallel)) {install.packages("doParallel")}
    if (!require(dplyr)) {install.packages("dplyr")}
    if (!require(future.apply)) {install.packages("future.apply")}
    
    # Required variables
    nrow_design <- nrow(design_matrix)
    name_res <- required_fx[1]
    name_t <- required_fx[2]
    seed <- as.integer(required_fx[3])
    
    if (parall == "Parallel") {
        # Detect
        if (missing(nclusters)) {
            ncluster <- detectCores() / 2
        } else {
            ncluster <- nclusters
        }
        # Create clusters and register them
        cl <- makeCluster(ncluster, type = "FORK")
        # Apply simulation
        clusterMap(cl, run_simulation, Row = 1:nrow_design,
                   MoreArgs = list(
                       name_results = name_res,
                       name_times = name_t,
                       design_matrix = design_matrix,
                       results_folder = folder,
                       seed = seed)
                   )
        # Stop parallelisation
        stopCluster(cl)
        
    } else if (parall == "forEach") {
        # Detect
        if (missing(nclusters)) {
            ncluster <- detectCores() / 2
        } else {
            ncluster <- nclusters
        }
        
        # Create clusters and register them
        cl <- makeCluster(ncluster, type = "FORK")
        registerDoParallel(cl)
        # Distribute rows into clusters
        rows_divided <- split(seq(nrow_design), 1:ncluster)
        # Export libraries, functions and variables
        # clusterExport(cl)
        #clusterExport(cl, required_fx, varlist = c("run_simulation", "SSD_mult_CRT"))
        # Parallelisation
        
        foreach(Rows = rows_divided, .errorhandling = "remove") %dopar% {
            run_simulation(Row = Rows, name_results = name_res, name_times = name_t,
                           design_matrix = design_matrix, results_folder = folder,
                           seed = seed)
        }
        # Stop clusters
        stopCluster(cl)
        stopImplicitCluster()
        
    } else if (parall == "future") {
        if (missing(nclusters)) {
            nclusters <- detectCores() / 2
        } else {
            nclusters <- nclusters
        }
        plan(multisession, workers = nclusters)
        rows_to_run <- 1:nrow_design
        future_lapply(rows_to_run, function(Row){
            run_simulation(Row, name_results = name_res, name_times = name_t,
                           design_matrix = design_matrix, results_folder = folder,
                           seed = seed)
        })
        
    }
    
}


# Collect results in a matrix---------------------------------------------
collect_results <- function(design_matrix, results_folder, finding, pair, name_results, test) {
    rows <-  seq(nrow(design_matrix))
    if (finding == "N2") {
        results_name <- ifelse(missing(name_results), "/ResultsN2Row", paste0("/", name_results))
        file_name <- "final_results_findN2"
    } else if (finding == "N1") {
        results_name <- "/ResultsN1Row"
        file_name <- "final_results_findN1"
    }
    
    if (test == "intersection-union") {
        new_matrix <- matrix(NA, ncol = 18, nrow = nrow(design_matrix))
        # extract results
        for (row_design in rows) {
            stored_result <- readRDS(paste0(results_folder, results_name, row_design, ".RDS"))
            n2 <- stored_result$n2
            n1 <- stored_result$n1
            median.BF12 <- median(stored_result[[11]]$results_H1[, "BF.12"])
            median.BF13 <- median(stored_result[[11]]$results_H1[, "BF.13"])
            median.BF14 <- median(stored_result[[11]]$results_H1[, "BF.14"])
            median.BF1c <- median(stored_result[[11]]$results_H1[, "BF.1c"])
            median.BF21 <- median(stored_result[[11]]$results_H2[, "BF.21"])
            median.BF2c <- median(stored_result[[11]]$results_H2[, "BF.2c"])
            median.BF31 <- median(stored_result[[11]]$results_H3[, "BF.31"])
            median.BF3c <- median(stored_result[[11]]$results_H3[, "BF.3c"])
            median.BF41 <- median(stored_result[[11]]$results_H4[, "BF.41"])
            median.BF4c <- median(stored_result[[11]]$results_H4[, "BF.4c"])
            eta.BF12 <- stored_result$Proportion.BF12
            eta.BF13 <- stored_result$Proportion.BF13
            eta.BF14 <- stored_result$Proportion.BF14
            eta.BF21 <- stored_result$Proportion.BF21
            eta.BF31 <- stored_result$Proportion.BF31
            eta.BF41 <- stored_result$Proportion.BF41
            new_matrix[row_design, ] <- c(median.BF12,
                                          median.BF13,
                                          median.BF14,
                                          median.BF1c, median.BF21, median.BF2c,
                                          median.BF31, median.BF3c, median.BF41,
                                          median.BF4c, eta.BF12, eta.BF13,
                                          eta.BF14, eta.BF21, eta.BF31, eta.BF41,
                                          n2, n1)
        }
        # create final data frame
        new_matrix <- as.data.frame(cbind(design_matrix, new_matrix))
        colnames(new_matrix) <- c(names(design_matrix), "median.BF12",
                                  "median.BF13", "median.BF14", "median.BF1c",
                                  "median.BF21",
                                  "median.BF2c", "median.BF31", "median.BF3c",
                                  "median.BF41", "median.BF4c", "eta.BF12",
                                  "eta.BF13", "eta.BF14", "eta.BF21", "eta.BF31",
                                  "eta.BF41", "n2.final", "n1.final")
        
    } else if (test == "homogeneity") {
        
    } else if (test == "omnibus") {
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
    
    # save results and return
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
