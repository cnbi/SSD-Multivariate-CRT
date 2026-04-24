########### SIMULATION FUNCTIONS ###############


# simulation function for second project
# name_results: characters. Name which will be use fo saving results
# name_times: characters. Name which will be use fo saving time results
run_simulation <- function(Row, name_results, name_times, design_matrix, results_folder){
    # Start time
    start_time <- Sys.time()
    
    # Actual simulation
    ssd_results <- SSD_mult_CRT(test = design_matrix[Row, "test"], 
                                effect_sizes = c(design_matrix[Row, "eff_size1"], design_matrix[Row, "eff_size2"]), 
                                n1 = design_matrix[Row, "n1"],
                                n2 = design_matrix[Row, "n2"], 
                                ndatasets = 500,
                                out_specific_ICC = c(design_matrix[Row, "out_specific_ICC"], 0.1), 
                                intersubj_between_outICC = design_matrix[Row, "intersubj_between_outICC"], 
                                intrasubj_between_outICC = design_matrix[Row, "intrasubj_between_outICC"],
                                pmp_thresh = design_matrix[Row, "pmp_thresh"], eta = design_matrix[Row, "eta"], 
                                fixed = as.character(design_matrix[Row, "fixed"]), max = 500, 
                                Bayes_pack = as.character(design_matrix[Row, "Bayes_pack"]),
                                master.seed = design_matrix[Row, "seed"],
                                difference = design_matrix[Row, "delta"]) # This is only for homogeneity
    
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
    return(NULL)
}

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
collect_results <- function(design_matrix, results_folder, finding, name_results, test,
                            file_name, rows, save) {
    if (missing(rows)) {
        rows <- seq(nrow(design_matrix))
    }
    if (finding == "N2") {
        results_name <- ifelse(missing(name_results), "/ResultsN2Row", paste0("/", name_results))
    } else if (finding == "N1") {
        results_name <- ifelse(missing(name_results), "/ResultsN1Row", paste0("/", name_results))
    }
    
    if (test == "omnibus") {
        new_matrix <- matrix(NA, ncol = 28, nrow = nrow(design_matrix))
        # extract results
        for (row_design in rows) {
            stored_result <- readRDS(paste0(results_folder, results_name, row_design, ".RDS"))
            n2 <- stored_result$n2
            n1 <- stored_result$n1
            median.BF14 <- median(stored_result$data$results_H1[, "BF.1c"])
            mean.PMP2.H1 <- mean(stored_result$data$results_H1[, "PMP.2"])
            mean.PMP3.H1 <- mean(stored_result$data$results_H1[, "PMP.3"])
            mean.PMP1.H1 <- mean(stored_result$data$results_H1[, "PMP.1"])
            median.BF24 <- median(stored_result$data$results_H2[, "BF.2c"])
            mean.PMP2.H2 <- mean(stored_result$data$results_H2[, "PMP.2"])
            mean.PMP3.H2 <- mean(stored_result$data$results_H2[, "PMP.3"])
            mean.PMP1.H2 <- mean(stored_result$data$results_H2[, "PMP.1"])
            median.BF34 <- median(stored_result$data$results_H3[, "BF.3c"])
            mean.PMP2.H3 <- mean(stored_result$data$results_H3[, "PMP.2"])
            mean.PMP3.H3 <- mean(stored_result$data$results_H3[, "PMP.3"])
            mean.PMP1.H3 <- mean(stored_result$data$results_H3[, "PMP.1"])
            median.BF41 <- median(stored_result$data$results_H4[, "BF.41"])
            median.BF42 <- median(stored_result$data$results_H4[, "BF.42"])
            median.BF43 <- median(stored_result$data$results_H4[, "BF.43"])
            mean.PMP2.H4 <- mean(stored_result$data$results_H4[, "PMP.2"])
            mean.PMP3.H4 <- mean(stored_result$data$results_H4[, "PMP.3"])
            mean.PMP1.H4 <- mean(stored_result$data$results_H4[, "PMP.1"])
            eta.PMP1 <- stored_result$Proportion.PMP1
            eta.PMP2 <- stored_result$Proportion.PMP2
            eta.PMP3 <- stored_result$Proportion.PMP3
            eta.PMP4 <- stored_result$Proportion.PMP4
            combined.PMP.H1 <- mean(stored_result$Combined.PMP.H1)
            combined.PMP.H2 <- mean(stored_result$Combined.PMP.H2)
            combined.PMP.H3 <- mean(stored_result$Combined.PMP.H3)
            PMP.H4 <- mean(stored_result$PMP.H4)
            new_matrix[row_design, ] <- c(median.BF14, mean.PMP2.H1, mean.PMP3.H1,
                                          mean.PMP1.H1, median.BF24, mean.PMP2.H2,
                                          mean.PMP3.H2, mean.PMP1.H2, median.BF34,
                                          mean.PMP2.H3, mean.PMP3.H3, mean.PMP1.H3,
                                          median.BF41, median.BF42, median.BF43,
                                          mean.PMP2.H4, mean.PMP3.H4, mean.PMP1.H4,
                                          eta.PMP1, eta.PMP2, eta.PMP3, eta.PMP4,
                                          combined.PMP.H1, combined.PMP.H2, combined.PMP.H3,
                                          PMP.H4,
                                          n2, n1)
        }
        # create final data frame
        new_matrix <- as.data.frame(cbind(design_matrix, new_matrix))
        colnames(new_matrix) <- c(names(design_matrix), "median.BF14", "mean.PMP2.H1",
                                  "mean.PMP3.H1", 
                                  "mean.PMP1.H1", "median.BF24", "mean.PMP2.H2", 
                                  "mean.PMP3.H2", "mean.PMP1.H2", "median.BF34", 
                                  "mean.PMP2.H3", "mean.PMP3.H3", "mean.PMP1.H3", 
                                  "median.BF41", "median.BF42", "median.BF43", 
                                  "mean.PMP2.H4", "mean.PMP3.H4", "mean.PMP1.H4", 
                                  "eta.PMP1", "eta.PMP2", "eta.PMP3", "eta.PMP4", 
                                  "combined.PMP.H1", "combined.PMP.H2", "combined.PMP.H3", 
                                  "PMP.H4", "n2.final", "n1.final")
        
    } else if (test == "homogeneity") {
        new_matrix <- matrix(NA, nrow = nrow(design_matrix), ncol = 10)
        # Extract results
        for (row_design in rows) {
            stored_result <- readRDS(paste0(results_folder, results_name, row_design, ".RDS"))
            n2 <- stored_result$n2
            n1 <- stored_result$n1
            median.BF1c <- median(stored_result$data$results_H1[, "BF.1c"])
            median.BFc1 <- median(stored_result$data$results_H1[, "BF.c1"])
            mean.PMP1 <- mean(stored_result$data$results_H1[, "PMP.1"])
            median.BF2c <- median(stored_result$data$results_H2[, "BF.2c"])
            median.BFc2 <- median(stored_result$data$results_H2[, "BF.c2"])
            mean.PMP2 <- mean(stored_result$data$results_H2[, "PMP.2"])
            eta.PMP1 <- stored_result$Proportion.PMP1
            eta.PMP2 <- stored_result$Proportion.PMP2
            new_matrix[row_design, ] <- c(median.BF1c, median.BFc1,
                                          mean.PMP1,
                                          median.BF2c, median.BFc2, mean.PMP2,
                                          eta.PMP1, eta.PMP2,
                                          n2, n1)
        }
        # create final data frame
        new_matrix <- as.data.frame(cbind(design_matrix, new_matrix))
        colnames(new_matrix) <- c(names(design_matrix), "median.BF1c", 
                                  "median.BFc1", "mean.PMP1", "median.BF2c",
                                  "median.BFc2", "mean.PMP2", "eta.PMP1",
                                  "eta.PMP2", "n2.final", "n1.final")
    } else if (test == "intersection-union") {
        new_matrix <- matrix(NA, ncol = 18, nrow = nrow(design_matrix))
        # extract results
        for (row_design in rows) {
            stored_result <- readRDS(paste0(results_folder, results_name, row_design, ".RDS"))
            n2 <- stored_result$n2
            n1 <- stored_result$n1
            median.BF12 <- median(stored_result$data$results_H1[, "BF.12"])
            median.BF13 <- median(stored_result$data$results_H1[, "BF.13"])
            median.BF14 <- median(stored_result$data$results_H1[, "BF.1c"])
            mean.PMP1 <- mean(stored_result$data$results_H1[, "PMP.1c"])
            median.BF21 <- median(stored_result$data$results_H2[, "BF.21"])
            median.BF2c <- median(stored_result$data$results_H2[, "BF.2c"])
            mean.PMP2 <- mean(stored_result$data$results_H2[, "PMP.2c"])
            median.BF31 <- median(stored_result$data$results_H3[, "BF.31"])
            median.BF3c <- median(stored_result$data$results_H3[, "BF.3c"])
            mean.PMP3 <- mean(stored_result$data$results_H3[, "PMP.3c"])
            median.BF41 <- median(stored_result$data$results_H4[, "BF.41"])
            mean.PMP4 <- mean(stored_result$data$results_H4[, "PMP.c"])
            eta.PMP1 <- stored_result$Proportion.PMP1
            eta.PMP2 <- stored_result$Proportion.PMP2
            eta.PMP3 <- stored_result$Proportion.PMP3
            eta.PMP4 <- stored_result$Proportion.PMP4
            new_matrix[row_design, ] <- c(median.BF12, median.BF13, median.BF14,
                                          mean.PMP1,
                                          median.BF21, median.BF2c, mean.PMP2, 
                                          median.BF31, median.BF3c, mean.PMP3,
                                          median.BF41, mean.PMP4, 
                                          eta.PMP1, eta.PMP2,
                                          eta.PMP3, eta.PMP4,
                                          n2, n1)
        }
        # create final data frame
        new_matrix <- as.data.frame(cbind(design_matrix, new_matrix))
        colnames(new_matrix) <- c(names(design_matrix), "median.BF14",
                                  "mean.PMP2.H1", "mean.PMP3.H1", "mean.PMP1",
                                  "median.BF21", "median.BF2c", "mean.PMP2",
                                  "median.BF31", "median.BF3c", "mean.PMP3",
                                  "median.BF41", "mean.PMP4", "eta.PMP1",
                                  "eta.PMP2", "eta.PMP3", "eta.PMP4",
                                  "n2.final", "n1.final")
        
    }
    
    # save results and return
    if (save == TRUE) {
        saveRDS(new_matrix, file = file.path(results_folder, paste0(file_name, ".RDS")))
    }
    
    return(new_matrix)
}

# Collect times in a matrix ----
collect_times <- function(design_matrix, results_folder, finding, name_results, test,
                          file_name, rows) {
    if (missing(rows)) {
        rows <- seq(nrow(design_matrix))
    }
    new_matrix <- matrix(NA, nrow = nrow(design_matrix), ncol = 1)
    if (finding == "N2") {
        results_name <- ifelse(missing(name_results), "/timeN2Row", paste0("/", name_results))
        file_name <- "final_times_findN2"
    } else if (finding == "N1") {
        results_name <- ifelse(missing(name_results), "/timeN1Row", paste0("/", name_results))
        file_name <- "final_times_findN1"
    }
    
    for (row_result in rows) {
        stored_result <- readRDS(paste0(results_folder, results_name, row_result, ".RDS"))
        new_matrix[row_result] <- stored_result
    }
    new_matrix <- as.data.frame(cbind(design_matrix, new_matrix))
    colnames(new_matrix) <- c(names(design_matrix), "total.time")
    saveRDS(new_matrix, file = file.path(results_folder, paste0(file_name, ".RDS")))
    
    return(new_matrix)
}

# Check simulation--------------------
library(tools)
missing_rows <- function(folder_path,
                         name_pattern = NULL,
                         check_numbers,
                         underscore = TRUE) {
    files_names <- list.files(folder_path)
    
    # Filter by name pattern
    if (!is.null(name_pattern)) {
        filtered_names <- files_names[grep(name_pattern, files_names)]
    }
    
    # Extract the number of row
    if (underscore) {
        row_numbers <- sapply(filtered_names, function(names){
            parts <- unlist(strsplit(names, "_"))
            number <- parts[length(parts)]
            
            # Remove extension
            number_only <- sub("\\.[^.]+$", "", number)
            return(number_only)
        })   
    } else if (underscore == FALSE) {
        
        row_numbers <- sapply(filtered_names, function(names){
            # Remove extension
            name_no_ext <- file_path_sans_ext(names)
            
            number <- regmatches(name_no_ext, gregexpr("\\d+$", name_no_ext))[[1]]
            
            return(number)
        })
    }
    
    
    difference <- setdiff(check_numbers, row_numbers)
    print(difference)
}

# Collect results for plots ==========================

collect_results_pl <- function(design_matrix, results_folder, finding, name_results, test,
                            file_name, rows, save, ndatasets) {
    if (missing(rows)) {
        rows <- seq(nrow(design_matrix))
    }
    if (finding == "N2") {
        results_name <- ifelse(missing(name_results), "/ResultsN2Row", paste0("/", name_results))
    } else if (finding == "N1") {
        results_name <- ifelse(missing(name_results), "/ResultsN1Row", paste0("/", name_results))
    }
    
    if (test == "omnibus") {
        new_matrix <- matrix(NA, ncol = 16, nrow = (nrow(design_matrix)*ndatasets))
        # extract results
        for (row_design in rows) {
            stored_result <- readRDS(paste0(results_folder, results_name, row_design, ".RDS"))
            n2 <- stored_result$n2
            n1 <- stored_result$n1
            BF1u <- stored_result$data$results_H1[, "BF.1u"]
            PMP1.H1 <- stored_result$data$results_H1[, "PMP.1"]
            BF2u <- stored_result$data$results_H2[, "BF.2u"]
            PMP2.H2 <- stored_result$data$results_H2[, "PMP.2"]
            BF3u <- stored_result$data$results_H3[, "BF.3u"]
            PMP3.H3 <- stored_result$data$results_H3[, "PMP.3"]
            BF4u <- stored_result$data$results_H4[, "BF.4u"]
            # combined.PMP.H1 <- mean(stored_result$Combined.PMP.H1)
            # combined.PMP.H2 <- mean(stored_result$Combined.PMP.H2)
            # combined.PMP.H3 <- mean(stored_result$Combined.PMP.H3)
            PMP.H4 <- stored_result$PMP.H4
            
            new_data <- cbind(BF1u, PMP1.H1, BF2u,
                              PMP2.H2, BF3u, PMP3.H3,
                              BF4u, PMP.H4,
                              n1, n2)
            new_matrix[(1 + (ndatasets*(row_design - 1))):(ndatasets*row_design), 7:16] <- new_data
        }
        # create final data frame
        new_design_matrix <- design_matrix[, c("eff_size1",  "eff_size2", 
                                               "out_specific_ICC", "intersubj_between_outICC",
                                               "intrasubj_between_outICC", "pmp_thresh")]
        
        new_matrix[, 1:6] <- new_design_matrix %>% slice(rep(1:n(), each = ndatasets)) %>% 
            as.matrix()
        colnames(new_matrix) <- c("eff_size1",  "eff_size2", 
                                  "out_specific_ICC", "intersubj_between_outICC",
                                  "intrasubj_between_outICC", "pmp_thresh",
                                  "BF.1u", "PMP.1",
                                  "BF.2u", 
                                  "PMP.2", "BF.3u", "PMP.3", 
                                  "BF.4u", "PMP.4", "n1.final", "n2.final")
        new_matrix <- as.data.frame(new_matrix) %>% 
            pivot_longer(
                cols = c( "BF.1u", "PMP.1", "BF.2u", "PMP.2", "BF.3u", "PMP.3", 
                          "BF.4u", "PMP.4"),
                names_to = c(".value", "hypothesis"), 
                names_pattern = "(BF|PMP)\\.(\\d)"
            )
        
    } else if (test == "homogeneity") {
        new_matrix <- matrix(NA, nrow = nrow(design_matrix)*ndatasets, ncol = 13)
        # Extract results
        for (row_design in rows) {
            stored_result <- readRDS(paste0(results_folder, results_name, row_design, ".RDS"))
            n2 <- stored_result$n2
            n1 <- stored_result$n1
            BF1c <- stored_result$data$results_H1[, "BF.1c"]
            BFc1 <- stored_result$data$results_H2[, "BF.2c"]
            PMP1 <- stored_result$data$results_H1[, "PMP.1"]
            PMP2 <- stored_result$data$results_H2[, "PMP.2"]
            new_data <- cbind(BF1c, BFc1,
                          PMP1, PMP2,
                          n1, n2)
            new_matrix[(1 + (ndatasets * (row_design - 1))):(ndatasets*row_design), 8:13] <- new_data
        }
        
        # create final data frame
        new_design_matrix <- design_matrix[, c("eff_size1",  "eff_size2", 
                                               "out_specific_ICC", "intersubj_between_outICC",
                                               "intrasubj_between_outICC", "pmp_thresh",
                                               "delta")]
        new_matrix[, 1:7] <- new_design_matrix %>% slice(rep(1:n(), each = ndatasets)) %>% as.matrix()
        colnames(new_matrix) <- c("eff_size1",  "eff_size2", 
                                   "out_specific_ICC", "intersubj_between_outICC",
                                   "intrasubj_between_outICC", "pmp_thresh",
                                   "delta", "BF1", "BF2",
                                          "PMP1", "PMP2", "n1.final", "n2.final")
        
        new_matrix <- as.data.frame(new_matrix) %>% 
            pivot_longer(
                cols = c("BF1", "BF2", "PMP1", "PMP2"),
                names_to = c(".value", "hypothesis"), 
                names_pattern = "(BF|PMP)(\\d)"
            )
        
    } else if (test == "intersection-union") {
        new_matrix <- matrix(NA, ncol = 16, nrow = nrow(design_matrix)*ndatasets)
        # extract results
        for (row_design in rows) {
            stored_result <- readRDS(paste0(results_folder, results_name, row_design, ".RDS"))
            n2 <- stored_result$n2
            n1 <- stored_result$n1
            BF1u <- stored_result$data$results_H1[, "BF.1u"]
            PMP1.H1 <- stored_result$data$results_H1[, "PMP.1c"]
            BF2u <- stored_result$data$results_H2[, "BF.2u"]
            PMP2.H2 <- stored_result$data$results_H2[, "PMP.2c"]
            BF3u <- stored_result$data$results_H3[, "BF.3u"]
            PMP3.H3 <- stored_result$data$results_H3[, "PMP.3c"]
            BF4u <- stored_result$data$results_H4[, "BF.4u"]
            PMP.H4 <- stored_result$data$results_H4[, "PMP.c"]
            new_matrix[(1 + (ndatasets*(row_design - 1))):(ndatasets*row_design), 7:16] <- cbind(BF1u, PMP1.H1, BF2u,
                                                                                             PMP2.H2, BF3u, PMP3.H3,
                                                                                             BF4u, PMP.H4,
                                                                                             n1, n2)
        }
        # create final data frame
        new_design_matrix <- design_matrix[, c("eff_size1",  "eff_size2", 
                                               "out_specific_ICC", "intersubj_between_outICC",
                                               "intrasubj_between_outICC", "pmp_thresh")]
        new_matrix[, 1:6] <- new_design_matrix %>% slice(rep(1:n(), each = ndatasets)) %>% 
            as.matrix()
        colnames(new_matrix) <- c("eff_size1",  "eff_size2", 
                                  "out_specific_ICC", "intersubj_between_outICC",
                                  "intrasubj_between_outICC", "pmp_thresh",
                                  "BF.1u", "PMP.1",
                                  "BF.2u", 
                                  "PMP.2", "BF.3u", "PMP.3", 
                                  "BF.4u", "PMP.4", "n1.final", "n2.final")
        new_matrix <- as.data.frame(new_matrix) %>% 
            pivot_longer(
                cols = c( "BF.1u", "PMP.1", "BF.2u", "PMP.2", "BF.3u", "PMP.3", 
                          "BF.4u", "PMP.4"),
                names_to = c(".value", "hypothesis"), 
                names_pattern = "(BF|PMP)\\.(\\d)"
            )
    }
    
    # save results and return
    if (save == TRUE) {
        saveRDS(new_matrix, file = file.path(results_folder, paste0(file_name, ".RDS")))
    }
    
    return(new_matrix)
}
