######################## MULTIVARIATE DATA GENERATION ##########################

#' ndatasets: Number of datasets to generate.
#' n1: Common cluster size
#' n2: Common number of clusters
#' effect_sizes: Vector with effect sizes or slopes.
#' out_specific_ICC: Outcome/endpoint-specific intraclass correlation coefficient.
#' intersubj_between_outICC: Intersubject between-endpoint/outcome intraclass correlation coefficient.
#' intrasubj_between_outICC: Intrasubject between-endpoint/outcomes intraclass correlation coefficient.

gen_multiv_data <- function(ndatasets, n1, n2, effect_sizes, out_specific_ICC, intersubj_between_outICC, 
                            intrasubj_between_outICC, n_outcomes, seed){
    
    
    # Matrix with intersubject intraclass correlation coefficients
    intersubj_iccs <- matrix(NA, n_outcomes, n_outcomes)
    diag_endpoint_icc <- c(out_specific_ICC, out_specific_ICC) #rho_0 of the two outcomes
    intersubj_iccs <- diag(diag_endpoint_icc)
    values_intersubj_mat <- intersubj_between_outICC #rho_1: order y2y1, y3y1, y3y2
    i <- 0
    for (column in 1:(n_outcomes - 1)) {
        for (row in 2:n_outcomes) {
            if (column != row) {
                i <- i + 1
                intersubj_iccs[row, column] <- values_intersubj_mat[i] #rho_1
                intersubj_iccs[column, row] <- values_intersubj_mat[i] #rho_1
            }
        }
    }
    y_names <- c("y1", "y2")
    colnames(intersubj_iccs) <- y_names
    rownames(intersubj_iccs) <- y_names
    
    # Matrix with intrasubject between-endpoints intraclass correlation coefficients
    intrasubj_icc <- matrix(NA, n_outcomes, n_outcomes)
    intrasubj_icc <- diag(1, n_outcomes)
    values_intrasubj_mat <- rep(intrasubj_between_outICC, (n_outcomes*(n_outcomes - 1)/2)) #rho_2: order y2y1, y3y1, y3y2
    i <- 0
    for (column in 1:(n_outcomes - 1)) {
        for (row in 2:n_outcomes) {
            if (column != row) {
                i <- i + 1
                intrasubj_icc[row, column] <- values_intrasubj_mat[i] #rho_2
                intrasubj_icc[column, row] <- values_intrasubj_mat[i] #rho_2
            }
        }
    }
    colnames(intrasubj_icc) <- y_names
    rownames(intrasubj_icc) <- y_names
    var_y <- rep(1, n_outcomes) # Marginal variance=total variance
    
    # Objects to save results
    output_multilevel <- vector(mode = "list", length = ndatasets)
    data_list <- vector(mode = "list", length = ndatasets)
    seeds <- vector(mode = "numeric", length = ndatasets)
    
    ## Treatment condition
    id_cluster <- rep(1:n2, each = n1)
    id_subj <- seq(1:(n1 * n2))
    condition <- rep(c(0, 1), each = n1 * n2 / 2)
    
    ## Covariance matrices specification 
    #Sigma e
    sigma_e <- diag(1 - diag(intersubj_iccs) * var_y)
    i <- 0
    for (column in 1:(n_outcomes - 1)) {
        for (row in 2:n_outcomes) {
            if (column != row) {
                i <- i + 1
                sigma_e[row, column] <- sqrt(var_y[row]) * sqrt(var_y[column]) * (intrasubj_icc[row, column] - intersubj_iccs[row, column])
                sigma_e[column, row] <- sqrt(var_y[row]) * sqrt(var_y[column]) * (intrasubj_icc[column, row] - intersubj_iccs[column, row])
            }
        }
    }
    
    ## Sigma u0
    sigma_u0 <- diag(diag(intersubj_iccs) * var_y)
    i <- 0
    for (column in 1:(n_outcomes - 1)) {
        for (row in 2:n_outcomes) {
            if (column != row) {
                i <- i + 1
                sigma_u0[row, column] <- sqrt(var_y[row]) * sqrt(var_y[column]) * intersubj_iccs[row, column]
                sigma_u0[column, row] <- sqrt(var_y[row]) * sqrt(var_y[column]) * intersubj_iccs[column, row]
            }
        }
    }
    sigma_e <- calibration_nonpos_def(sigma_e) # to avoid non-positive definite covariance matrix
    sigma_u0 <- calibration_nonpos_def(sigma_u0) # to avoid non-positive definite covariance matrix
    
    # Generate new seed in case of singular convergence
    max_attempts <- 10
    attempt <- 1
    success <- FALSE

    # Random effects
    while (attempt <= max_attempts && success == FALSE) {
        try({
            for (iteration in seq(ndatasets)) {
                seeds[iteration] <- (iteration + seed) * iteration * attempt
                print(seeds[iteration])
                set.seed(seeds[iteration])
                e <- MASS::mvrnorm(n1*n2, rep(0, n_outcomes), sigma_e, tol = 1e-04)
                u0 <- MASS::mvrnorm(n2, rep(0, n_outcomes), sigma_u0, tol = 1e-04)
                u0 <- do.call(rbind, replicate(n1, u0, simplify = FALSE))
                
                z_bar <- mean(condition)
                if (n_outcomes == 2) {
                    y1 <- 0 + effect_sizes[1] * (condition - z_bar) + u0[, 1] + e[, 1]
                    y2 <- 0 + effect_sizes[2] * (condition - z_bar) + u0[, 2] + e[, 2]
                    my_data <- cbind(id_subj, id_cluster, y1, y2, condition)
                } else if (n_outcomes == 3) {
                    y1 <- 0 + effect_sizes[1] * (condition - z_bar) + u0[, 1] + e[, 1]
                    y2 <- 0 + effect_sizes[2] * (condition - z_bar) + u0[, 2] + e[, 2]
                    y3 <- 0 + effect_sizes[3] * (condition - z_bar) + u0[, 3] + e[, 3]
                    my_data <- cbind(id_subj, id_cluster, y1, y2, y3, condition)
                }
                colnames(my_data)[2] <- "cluster"
                data_list[[iteration]] <- as.data.frame(my_data)
            }
            # EM
            print("EM begins")
            if (n_outcomes == 2) {
                estimations <- Map(EM.estim2, data_list, list(as.formula('y1 ~ condition + (1|cluster)')), list(as.formula('y2 ~ condition + (1|cluster)')), list(500))
            } else if (n_outcomes == 3) {
                estimations <- Map(EM.estim3, data_list, list(as.formula('y1 ~ condition + (1|cluster)')), list(as.formula('y2 ~ condition + (1|cluster)')), list(as.formula('y3 ~ condition + (1|cluster)')), list(500))
            }
            success <- TRUE
        }, silent = TRUE)
        if (success == FALSE) {
            warning("Attempt", attempt, "failed. Regenerating data")
            attempt <- attempt + 1
        }
    }
    if (success == FALSE) {
        
        stop("Data generation failed after reaching maximum attempts")
        
    }
    
    # estimations <- EM.estim2(as.data.frame(my_data), as.formula('y1 ~ condition'), as.formula('y2 ~ condition'), maxiter = 500, verbose = TRUE)
    
    # Extract the elements for outcomes
    fixed_eff <- Map(extract_fix_eff, estimations, list(n_outcomes))
    
    # Extract variance-covariance matrix of parameters of interest
    var_cov <- Map(extract_var_cov, estimations, list(n_outcomes))
    
    # Calculate ICCs
    ICCs <- Map(calc_ICCs, estimations, list(n_outcomes))
    print("Data generation done!")
    return(list(estimations = fixed_eff,
                Sigma = var_cov,
                ICCs = ICCs,
                seeds = seeds))
}
