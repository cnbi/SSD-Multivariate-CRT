######################## MULTIVARIATE DATA GENERATION ##########################

#' ndatasets: Number of datasets to generate.
#' n1: Common cluster size
#' n2: Total number of clusters
#' effect_sizes: Vector with effect sizes or slopes.
#' out_specific_ICCs: Vector with outcome/endpoint-specific intraclass correlation coefficient.
#' intersubj_between_outICC: Intersubject between-endpoint/outcome intraclass correlation coefficient.
#' intrasubj_between_outICC: Intrasubject between-endpoint/outcomes intraclass correlation coefficient.

# TODO: The following code only works for n_outcomes =  2

gen_multiv_data <- function(ndatasets, n1, n2, effect_sizes, out_specific_ICCs, intersubj_between_outICC, 
                            intrasubj_between_outICC, n_outcomes, master.seed, homogeneity = FALSE){
    
    if (homogeneity == TRUE) {
        marginal_variances <- c(30, 30)
    } else{
        set.seed(master.seed)
        marginal_variances <- runif(n_outcomes, 10, 50)
    }
    
    # Calculate the variance components
    var_u0_scaled <- out_specific_ICCs * marginal_variances
    var_e_scaled <- marginal_variances * (1 - out_specific_ICCs)
    
    # Defining variance-covariance matrices of random effects
    sigma_u0 <- diag(var_u0_scaled)
    sigma_e <- diag(var_e_scaled)
    
    # Covariances components
    cov_e <- matrix(0, n_outcomes, n_outcomes)
    for (column in 1:(n_outcomes - 1)) {
        for (row in 2:n_outcomes) {
            if (column != row) {
                # At cluster level
                sigma_u0[row, column] <- intersubj_between_outICC * sqrt(marginal_variances[column] * marginal_variances[row])#rho_1
                sigma_u0[column, row] <- intersubj_between_outICC * sqrt(marginal_variances[column] * marginal_variances[row]) #rho_1
                
                total_covariance <- intrasubj_between_outICC * sqrt(marginal_variances[row] * marginal_variances[column])
                
                # At individual level
                sigma_e[row, column] <- total_covariance - sigma_u0[row, column]
                sigma_e[column, row] <- total_covariance - sigma_u0[column, row]
            }
        }
    }
    y_names <- c("y1", "y2")
    colnames(sigma_u0) <- y_names
    rownames(sigma_u0) <- y_names
    colnames(sigma_e) <- y_names
    rownames(sigma_e) <- y_names
    
    # ID
    cluster <- rep(1:n2, each = n1)
    id_subj <- seq(1:(n1 * n2))
    
    # Treatment condition
    condition_cluster <- rep(c(0, 1), each = n2 / 2) # randomisation of clusters
    condition <- rep(condition_cluster, each = n1)
    
    # Scaled effect sizes
    scaled_effects <- effect_sizes * sqrt(marginal_variances)
    control_means <- 0
    
    # Objects to save results
    output_multilevel <- vector(mode = "list", length = ndatasets)
    data_list <- vector(mode = "list", length = ndatasets)
    seeds <- seq(ndatasets) + master.seed
    # Random effects
    for (iteration in seq(ndatasets)) {
        set.seed(seeds[iteration])
        e <- MASS::mvrnorm(n1 * n2, rep(0, n_outcomes), sigma_e)
        u0 <- MASS::mvrnorm(n2, rep(0, n_outcomes), sigma_u0)
        u0_y1 <- u0[, 1]
        u0_y2 <- u0[, 2]
        u0_y1 <- rep.int(u0_y1,times = rep(n1, n2))
        u0_y2 <- rep.int(u0_y2,times = rep(n1, n2))
        e_y1 <- e[, 1]
        e_y2 <- e[, 2]
        y1 <- control_means[1] + scaled_effects[1] * condition + u0_y1 + e_y1
        y2 <- control_means[2] + scaled_effects[2] * condition + u0_y2 + e_y2
        my_data <- cbind(id_subj, cluster, condition, y1, y2)
        data_list[[iteration]] <- as.data.frame(my_data)
    }
    
    # Multilevel SEM
    model <- "
    level: 1
    y1 ~~ y2
    
    level: 2
    y1 + y2 ~ condition
    y1 ~~ y2
    "
    
    # output_multilevel <- semList(model = model, dataList = data_list, cluster = "cluster", 
    #                       store.failed = TRUE, iseed = seed, show.progress = TRUE)
    
    output_multilevel <- Map(lavaan::sem, model = list(model), data = data_list, cluster = "cluster")
    
    #I can use clusterApply for a parallel version
    
    # Extract info
    fixed_eff <- Map(extract_fix_eff, output_multilevel, list(n_outcomes))
    var_cov <- Map(extract_var_cov, output_multilevel, list(n_outcomes))
    # Calculate ICCs
    ICCs <- Map(calc_ICCs, output_multilevel, list(n_outcomes))

    print("Data generation done!")
    #return(output_multilevel)
    return(list(estimations = fixed_eff, #Vector with unstd. and std intercepts and beta1
                Sigma = var_cov, #Matrix
                ICCs = ICCs, #List with matrices and a vector containing the ICCs
                seeds = seeds)) #Vector with seeds
}


# Version 2---------------------------------------------------------------------
# gen_multiv_data2 <- function(ndatasets, n1, n2, effect_sizes, out_specific_ICCs, intersubj_between_outICC, 
#                             intrasubj_between_outICC, n_outcomes, seed){
#     
#     # Matrix with intersubject intraclass correlation coefficients
#     intersubj_iccs <- matrix(NA, n_outcomes, n_outcomes)
#     diag_endpoint_icc <- c(out_specific_ICCs, out_specific_ICCs) #rho_0 of the two outcomes
#     intersubj_iccs <- diag(diag_endpoint_icc)
#     values_intersubj_mat <- intersubj_between_outICC #rho_1: order y2y1, y3y1, y3y2
#     i <- 0
#     for (column in 1:(n_outcomes - 1)) {
#         for (row in 2:n_outcomes) {
#             if (column != row) {
#                 i <- i + 1
#                 intersubj_iccs[row, column] <- values_intersubj_mat[i] #rho_1
#                 intersubj_iccs[column, row] <- values_intersubj_mat[i] #rho_1
#             }
#         }
#     }
#     y_names <- c("y1", "y2")
#     colnames(intersubj_iccs) <- y_names
#     rownames(intersubj_iccs) <- y_names
#     
#     # Matrix with intrasubject between-endpoints intraclass correlation coefficients
#     intrasubj_icc <- matrix(NA, n_outcomes, n_outcomes)
#     intrasubj_icc <- diag(1, n_outcomes)
#     values_intrasubj_mat <- rep(intrasubj_between_outICC, (n_outcomes*(n_outcomes - 1)/2)) #rho_2: order y2y1, y3y1, y3y2
#     i <- 0
#     for (column in 1:(n_outcomes - 1)) {
#         for (row in 2:n_outcomes) {
#             if (column != row) {
#                 i <- i + 1
#                 intrasubj_icc[row, column] <- values_intrasubj_mat[i] #rho_2
#                 intrasubj_icc[column, row] <- values_intrasubj_mat[i] #rho_2
#             }
#         }
#     }
#     colnames(intrasubj_icc) <- y_names
#     rownames(intrasubj_icc) <- y_names
#     var_y <- rep(1, n_outcomes) # Marginal variance=total variance
#     
#     # Objects to save results
#     output_multilevel <- vector(mode = "list", length = ndatasets)
#     data_list <- vector(mode = "list", length = ndatasets)
#     seeds <- vector(mode = "numeric", length = ndatasets)
#     
#     ## Treatment condition
#     cluster <- rep(1:n2, each = n1)
#     id_subj <- seq(1:(n1 * n2))
#     condition <- rep(c(0, 1), each = n1 * n2 / 2)
#     
#     ## Covariance matrices specification 
#     #Sigma e
#     sigma_e <- diag(1 - diag(intersubj_iccs) * var_y)
#     i <- 0
#     for (column in 1:(n_outcomes - 1)) {
#         for (row in 2:n_outcomes) {
#             if (column != row) {
#                 i <- i + 1
#                 sigma_e[row, column] <- sqrt(var_y[row]) * sqrt(var_y[column]) * (intrasubj_icc[row, column] - intersubj_iccs[row, column])
#                 sigma_e[column, row] <- sqrt(var_y[row]) * sqrt(var_y[column]) * (intrasubj_icc[column, row] - intersubj_iccs[column, row])
#             }
#         }
#     }
#     
#     ## Sigma u0
#     sigma_u0 <- diag(diag(intersubj_iccs) * var_y)
#     i <- 0
#     for (column in 1:(n_outcomes - 1)) {
#         for (row in 2:n_outcomes) {
#             if (column != row) {
#                 i <- i + 1
#                 sigma_u0[row, column] <- sqrt(var_y[row]) * sqrt(var_y[column]) * intersubj_iccs[row, column]
#                 sigma_u0[column, row] <- sqrt(var_y[row]) * sqrt(var_y[column]) * intersubj_iccs[column, row]
#             }
#         }
#     }
#     sigma_e <- calibration_nonpos_def(sigma_e) # to avoid non-positive definite covariance matrix
#     sigma_u0 <- calibration_nonpos_def(sigma_u0) # to avoid non-positive definite covariance matrix
#     
#     
#     # Variance-covariance matrix
#     var_cov_matrix <- var_cov(n2 = n2, n1 = n1, Sigma_e = sigma_e, Sigma_u = sigma_u0)
#     # var_cov_matrix <- calibration_nonpos_def(var_cov_matrix) # to avoid non-positive definite covariance matrix
#     print("Variance-covariance matrix of treatment effects (betas):")
#     print(var_cov_matrix)
#     
#     # Generate bivariate normal samples of treatment effects
#     # Generate new seed in case of singular convergence
#     set.seed(seed)  # for reproducibility
#     sampled_effects <- MASS::mvrnorm(n = 5000, mu = effect_sizes, Sigma = var_cov_matrix)
#     
#     print("First samples of treatment effects from bivariate normal distribution:")
#     print(head(sampled_effects))
#     
#     print("Data generation done!")
#     return(list(estimations = sampled_effects,
#                 Sigma = var_cov_matrix,
#                 seeds = seed))
# }
# 
# # Version 3---------------------------------------------------------------------
# gen_multiv_data3 <- function(ndatasets, n1, n2, effect_sizes, out_specific_ICCs, intersubj_between_outICC, 
#                              intrasubj_between_outICC, n_outcomes, seed){
#     
#     # Matrix with intersubject intraclass correlation coefficients
#     intersubj_iccs <- matrix(NA, n_outcomes, n_outcomes)
#     diag_endpoint_icc <- c(out_specific_ICCs, out_specific_ICCs) #rho_0 of the two outcomes
#     intersubj_iccs <- diag(diag_endpoint_icc)
#     values_intersubj_mat <- intersubj_between_outICC #rho_1: order y2y1, y3y1, y3y2
#     i <- 0
#     for (column in 1:(n_outcomes - 1)) {
#         for (row in 2:n_outcomes) {
#             if (column != row) {
#                 i <- i + 1
#                 intersubj_iccs[row, column] <- values_intersubj_mat[i] #rho_1
#                 intersubj_iccs[column, row] <- values_intersubj_mat[i] #rho_1
#             }
#         }
#     }
#     y_names <- c("y1", "y2")
#     colnames(intersubj_iccs) <- y_names
#     rownames(intersubj_iccs) <- y_names
#     
#     # Matrix with intrasubject between-endpoints intraclass correlation coefficients
#     intrasubj_icc <- matrix(NA, n_outcomes, n_outcomes)
#     intrasubj_icc <- diag(1, n_outcomes)
#     values_intrasubj_mat <- rep(intrasubj_between_outICC, (n_outcomes*(n_outcomes - 1)/2)) #rho_2: order y2y1, y3y1, y3y2
#     i <- 0
#     for (column in 1:(n_outcomes - 1)) {
#         for (row in 2:n_outcomes) {
#             if (column != row) {
#                 i <- i + 1
#                 intrasubj_icc[row, column] <- values_intrasubj_mat[i] #rho_2
#                 intrasubj_icc[column, row] <- values_intrasubj_mat[i] #rho_2
#             }
#         }
#     }
#     colnames(intrasubj_icc) <- y_names
#     rownames(intrasubj_icc) <- y_names
#     var_y <- rep(1, n_outcomes) # Marginal variance=total variance
#     
#     # Objects to save results
#     output_multilevel <- vector(mode = "list", length = ndatasets)
#     data_list <- vector(mode = "list", length = ndatasets)
#     seeds <- vector(mode = "numeric", length = ndatasets)
#     
#     ## Treatment condition
#     id_cluster <- rep(1:n2, each = n1)
#     id_subj <- seq(1:(n1 * n2))
#     condition <- rep(c(0, 1), each = n1 * n2 / 2)
#     
#     ## Covariance matrices specification 
#     #Sigma e
#     sigma_e <- diag(1 - diag(intersubj_iccs) * var_y)
#     i <- 0
#     for (column in 1:(n_outcomes - 1)) {
#         for (row in 2:n_outcomes) {
#             if (column != row) {
#                 i <- i + 1
#                 sigma_e[row, column] <- sqrt(var_y[row]) * sqrt(var_y[column]) * (intrasubj_icc[row, column] - intersubj_iccs[row, column])
#                 sigma_e[column, row] <- sqrt(var_y[row]) * sqrt(var_y[column]) * (intrasubj_icc[column, row] - intersubj_iccs[column, row])
#             }
#         }
#     }
#     
#     ## Sigma u0
#     sigma_u0 <- diag(diag(intersubj_iccs) * var_y)
#     i <- 0
#     for (column in 1:(n_outcomes - 1)) {
#         for (row in 2:n_outcomes) {
#             if (column != row) {
#                 i <- i + 1
#                 sigma_u0[row, column] <- sqrt(var_y[row]) * sqrt(var_y[column]) * intersubj_iccs[row, column]
#                 sigma_u0[column, row] <- sqrt(var_y[row]) * sqrt(var_y[column]) * intersubj_iccs[column, row]
#             }
#         }
#     }
#     sigma_e <- calibration_nonpos_def(sigma_e) # to avoid non-positive definite covariance matrix
#     sigma_u0 <- calibration_nonpos_def(sigma_u0) # to avoid non-positive definite covariance matrix
#     
#     
#     # Variance-covariance matrix
#     var_cov_matrix <- var_cov2(n2 = n2, n1 = n1, Sigma_e = sigma_e, Sigma_u = sigma_u0, z_bar = 0.5)
#     # var_cov_matrix <- calibration_nonpos_def(var_cov_matrix) # to avoid non-positive definite covariance matrix
#     print("Variance-covariance matrix of treatment effects (betas):")
#     print(var_cov_matrix)
#     
#     # Generate bivariate normal samples of treatment effects
#     # Generate new seed in case of singular convergence
#     set.seed(seed)  # for reproducibility
#     sampled_effects <- MASS::mvrnorm(n = 5000, mu = effect_sizes, Sigma = var_cov_matrix)
#     
#     print("First samples of treatment effects from bivariate normal distribution:")
#     print(head(sampled_effects))
#     
#     print("Data generation done!")
#     return(list(estimations = sampled_effects,
#                 Sigma = var_cov_matrix,
#                 seeds = seed))
# }