######################## HELPERS #################

# Calculate effective sample size -------------------------------
# 
# Source: Hox, J. J., Moerbeek, M., & van de Schoot, R. (2017). 
## Multilevel Analysis: Techniques and Applications (3rd ed.). Routledge. https://doi.org/10.4324/9781315650982
effective_sample <- function(n1, n2, data, n_outcomes){
    icc <- data$emp_rho0
    effective_n <- (n1*n2) / (1 + (n1 - 1) * icc)
    return(effective_n)
}


# Negative log-likelihood ------------------------
#
neg_loglik <- function(theta, n2, Y, X, ID, n1s, outcomes) {
    neg_loglik <- -loglik(theta, n2, Y, X, ID, n1s, outcomes)
    return(neg_loglik)
}


# Extract results ----------------
#
extract_res <- function(x, number) {
    #results <- ifelse(is.null(x[[number]]), NaN, x[[number]])
    results <- x[[number]]
    return(results)
}


# Rounding half away from zero -----------------------------
#
# Source: https://stackoverflow.com/questions/66600344/commercial-rounding-in-r-i-e-always-round-up-from-5/66600470#66600470
round2 <- function(number, decimals = 0) {
    sign_number <- sign(number)
    number <- abs(number) * 10^decimals
    number <- number + 0.5 + sqrt(.Machine$double.eps)
    number <- trunc(number)
    number <- number / 10 ^ decimals
    number * sign_number
}


# Extract effect sizes -----------------------------------
#
extract_fix_eff <- function(data, n_outcomes){
    if (n_outcomes == 2) {
        # Extract treatment effect
        est_beta1 <- lavaan::lavInspect(data, "est")$cluster$beta[1:2, 3]
        st_beta1 <- lavaan::lavInspect(data, "std")$cluster$beta[1:2, 3]
        names(est_beta1) <- c("y1", "y2")
        names(st_beta1) <- c("y1", "y2")
        
        # Extract intercepts
        est_intercepts <- lavaan::lavInspect(data, "est")$cluster$alpha[1:2, ]
        st_intercepts <- lavaan::lavInspect(data, "std")$cluster$alpha[1:2, ]
        names(est_intercepts) <- c("y1", "y2")
        names(st_intercepts) <- c("y1", "y2")
        
        fixed_eff <- c(est_intercepts, est_beta1, st_intercepts, st_beta1)
        names(fixed_eff) <- c("Intercept1", "Intercept2", "Treatment1", "Treatment2",
                              "st.Intercept1", "st.Intercept2", "st.Treatment1", "st.Treatment2")
        
    } else if (n_outcomes == 3) {
        #TODO: Change code to extract for three outcomes
        # Extract the elements for outcome 1
        estimates_y1 <- data$theta$zeta[1:2]
        
        # Extract the elements for outcome 2
        estimates_y2 <- data$theta$zeta[3:4]
        
        # Extract the elements for outcome 3
        estimates_y3 <- data$theta$zeta[5:6]
        
        fixed_eff <- c(estimates_y1, estimates_y2, estimates_y3)
        names(fixed_eff) <- c("Intercept1", "Slope1", "Intercept2", "Slope2",
                              "Intercept3", "Slope3")
    }
    return(fixed_eff)
}


# Extract variance-covariance matrix --------------------------------------------
#
extract_var_cov <- function(lavaan_output, n_outcomes) {
    if (n_outcomes == 2) {
        var_cov_matrix <- lavaan::lavInspect(lavaan_output, "vcov")[4:5, 4:5]
    } else if (n_outcomes == 3) {
        param_interest <- c("y1~condition.l2", "y2~condition.l2", "y3~condition.l2")
        var_cov_matrix <- lavaan::lavInspect(lavaan_output, "vcov")[param_interest, param_interest]
    }
    return(var_cov_matrix)
}


# Calculate empirical ICCs ---------------------------------
# lavaan object
calc_ICCs <- function(lavaan_output, n_outcomes) {
    emp_rho0 <- vector(mode = "numeric", length = n_outcomes)  # Empirical Outcome-specific ICC
    emp_rho1 <- matrix(NA, n_outcomes, n_outcomes)  # Empirical intersubject between-outcome ICC
    emp_rho2 <- matrix(NA, n_outcomes, n_outcomes)  # Empirical intrasubject between-outcome ICC
    
    # Extract variance-covariance of random effects
    variance_cov_list <- lavaan::lavInspect(lavaan_output, "cov.ov")
    
    # Calculate empirical ICCs
    emp_Sigma_u0 <- variance_cov_list$cluster
    emp_Sigma_e <- variance_cov_list$within
    
    # Outcome-specific ICC
    for (y in 1:n_outcomes) {
        emp_rho0[y] <- emp_Sigma_u0[y, y]/(emp_Sigma_u0[y, y] + emp_Sigma_e[y, y])
    }
    if (n_outcomes == 2) {
        names(emp_rho0) <- c("y1", "y2")
    } else if (n_outcomes == 3) {
        names(emp_rho0) <- c("y1", "y2", "y3")
    }
    
    # Intersubject between-outcome ICC and Intrasubject between-outcome ICC
    for (y in 1:(n_outcomes - 1)) {
        for (y_prime in (y + 1):n_outcomes) {
            denominator <- sqrt((emp_Sigma_u0[y, y] + emp_Sigma_e[y, y]) * (emp_Sigma_u0[y_prime, y_prime] + emp_Sigma_e[y_prime, y_prime]))
            emp_rho1[y, y_prime] <- emp_Sigma_u0[y, y_prime]/denominator
            emp_rho1[y_prime, y] <- emp_rho1[y, y_prime]
            emp_rho2[y, y_prime] <- (emp_Sigma_u0[y, y_prime] + emp_Sigma_e[y, y_prime]) / denominator
            emp_rho2[y_prime, y] <- emp_rho2[y, y_prime]
        }
    }
    
    return(list(emp_rho0 = emp_rho0,
                emp_rho1 = emp_rho1,
                emp_rho2 = emp_rho2))
}

# Create the hypothesis to test------------------------------------
#
hypothesis_maker <- function(names_h, difference, constrain) {
    #estimates <- sort(estimates, decreasing = TRUE)
    hypothesis <- paste0(names_h[1],"-", names_h[2], constrain, as.character(difference))
    return(hypothesis)
}


#Calibration method for non-positive definite covariance matrix in multivariate data analysis---------------
# Source of the method: Chao Huang, Daniel Farewell, Jianxin Pan, A calibration method for non-positive definite covariance 
# matrix in multivariate data analysis, Journal of Multivariate Analysis, Volume 157,2017,Pages 45-52,ISSN 0047-259X,
# https://doi.org/10.1016/j.jmva.2017.03.001.
calibration_nonpos_def <- function(cov_matrix) {
    # Perform eigenvalue decomposition
    eigen_decomp <- eigen(cov_matrix)
    
    # Extract eigenvalues and eigenvectors
    eigenvalues <- eigen_decomp$values
    eigenvectors <- eigen_decomp$vectors
    
    # Check if any eigenvalue is non-positive and adjust the eigenvalues if necessary
    if (any(eigenvalues <= 0)) {
        # Adjust non-positive eigenvalues by adding a small constant to make them positive
        eigenvalues <- pmax(eigenvalues, 1e-4)  # Add a small constant (1e-6) to non-positive eigenvalues
        # Reconstruct the covariance matrix using the adjusted eigenvalues
        corrected_matrix <- eigenvectors %*% diag(eigenvalues) %*% t(eigenvectors)
        return(corrected_matrix)
    } else {
        return(cov_matrix)  # Return original matrix if it's already positive definite
    }
}

# Print hypotheses ---------------------------------------------
print_hypotheses <- function(list_hypo) {
    cat("Hypotheses:", "\n")
    for (hypoth in 1:length(list_hypo)) {
        cat(paste0("    H", hypoth, ":"), list_hypo[[hypoth]], "\n")
    }
}

# Compute variance-covariance matrix--------------------------------------------
# From Yang, et al. (2022, p. 1296). Estimation of treatment effects is carried out
# with the feasible generalised least square, whose large sample variance is given 
# by Equation 3 in the mentioned paper.

var_cov <- function(n2, n1, Sigma_e, Sigma_u) {
    P <- ncol(Sigma_e)
    Sigma_e_inv <- solve(Sigma_e)
    sum_Sigmas_inv <- solve(Sigma_e + n1 * Sigma_u)
    I_n1 <- diag(n1)                                # n1,n1 identity matrix
    J_n1 <- matrix(1, n1, n1)                       # n1,n1 matrix of ones
    
    V_inv_second <- (1/n1) * sum_Sigmas_inv - Sigma_e_inv
    V_i_inv <- kronecker(I_n1, Sigma_e_inv) + kronecker(J_n1, V_inv_second)
    
    # Design matrices for each treatment condition
    W_treatment <-  kronecker(rep(1, n1), cbind(diag(P), diag(P) * 1)) # [I_K, I_K]
    w_control <- kronecker(rep(1, n1), cbind(diag(P), diag(P) * 0)) # [I_K, 0]
    
    # Information matrices by treatment condition and total
    U_treatment <- t(W_treatment) %*% V_i_inv %*% W_treatment
    U_control <- t(w_control) %*% V_i_inv %*% w_control
    U <- (U_treatment * (n2/2)) + (U_control * (n2/2)) # Assuming equal allocation of clusters into treatment conditions
    cat("Before calibration")
    print(U)
    U <- calibration_nonpos_def2(U)
    cat("After calibration")
    print(U)
    
    var_parameters <- solve(U)
    
    # Extract variances and covariances
    var_cov_beta <- var_parameters[(P + 1):(2 * P), (P + 1):(2 * P)]
    return(var_cov_beta)
}

# Calibrate non-positive definite matrix 2
# Only change the negative eigenvalues
calibration_nonpos_def2 <- function(cov_matrix) {
    eigen_decomp <- eigen(cov_matrix)
    eigenvalues <- eigen_decomp$values
    eigenvectors <- eigen_decomp$vectors
    
    # Replace only non-positive eigenvalues with a small positive constant
    eigenvalues[eigenvalues <= 0] <- 1e-4
    
    corrected_matrix <- eigenvectors %*% diag(eigenvalues) %*% t(eigenvectors)
    return(corrected_matrix)
}

# Compute variance-covariance matrix2--------------------------------------------
# From Yang, et al. (2022, p. 1296). Estimation of treatment effects is carried out
# with the feasible generalised least square, whose large sample variance is given 
# by Equation 3 in the mentioned paper.

var_cov2 <- function(n2, n1, Sigma_e, Sigma_u, z_bar) {
    P <- ncol(Sigma_e) # Number of outcomes
    Sigma_e_inv <- solve(Sigma_e)
    sum_Sigmas_inv <- solve_matrix(Sigma_e + n1 * Sigma_u)
    I_n1 <- diag(n1)                                # n1,n1 identity matrix
    J_n1 <- matrix(1, n1, n1)                       # n1,n1 matrix of ones
    
    V_inv_second_term <- (1/n1) * sum_Sigmas_inv - Sigma_e_inv
    V_i_inv <- kronecker(I_n1, Sigma_e_inv) + kronecker(J_n1, V_inv_second_term)
    
    # Design matrices for each treatment condition
    W_treatment <-  kronecker(rep(1, n1), cbind(diag(P), diag(P) * (1 - z_bar)))
    w_control <- kronecker(rep(1, n1), cbind(diag(P), diag(P) * (0 - z_bar)))
    
    # Information matrices by treatment condition and total
    U_treatment <- t(W_treatment) %*% V_i_inv %*% W_treatment
    U_control <- t(w_control) %*% V_i_inv %*% w_control
    U <- (U_treatment * (n2/2)) + (U_control * (n2/2)) # Assuming equal allocation of clusters into treatment conditions
    
    cat("Before calibration")
    print(U)
    
    # Ensuring positive definiteness
    U_eigen_values <- eigen(U)$values
    if(any(U_eigen_values <= 0)) {
        U <- U + diag(2*P) * (abs(min(U_eigen_values)) + 1e-10)
    }
    
    cat("After calibration")
    print(U)
    
    var_parameters <- solve(U)
    
    # Extract variances and covariances
    var_cov_beta <- var_parameters[(P + 1):(2 * P), (P + 1):(2 * P)]
    return(var_cov_beta)
}

# Matrix inversion with ridge regularisation-----------------------------------
solve_matrix <- function(Matrix, tol = 1e-12) {
    cond_number <- kappa(Matrix, exact = TRUE)
    if (cond_number > 1/tol) {
        lambda <- max(tol, 1/cond_number)
        Matrix_reg <- Matrix + diag(nrow(Matrix)) * lambda
        return(solve(Matrix_reg))
    } else {
        return(solve(Matrix))
    }
}


# Input validation -----------------------------------------------------------.

check_inputs <- function(test, effect_sizes, n1 = 15, n2 = 30, ndatasets = 1000, out_specific_ICC, 
                         intersubj_between_outICC, intrasubj_between_outICC,
                         pmp_thresh = 0.9, eta = 0.8, fixed = "n1", difference = 0.2, max,
                         master.seed, Bayes_pack){
    
    numeric_input <-  c(effect_sizes, n1, n2, ndatasets, out_specific_ICC, 
                        intersubj_between_outICC, intrasubj_between_outICC, pmp_thresh, 
                        eta, max)
    
    if (!all(sapply(numeric_input, is.numeric))) {
        stop("All arguments, except 'fixed' and 'Bayes_pack', must be numeric")
    }
    
    if (n2 %% 2 > 0) stop("The number of clusters must be even")
    if (eta > 1 || eta < 0) {
        stop("The probability of exceeding the threshold must be between 0 and 1")
    }
    if (!is.character(fixed) || !fixed %in% c("n1", "n2")) stop("Fixed can only be a character indicating n1 or n2.")
    if ((test == "homogeneity") && (effect_sizes[1] < effect_sizes[2]))
        stop("Effect size 1 must be larger than effect size 2")
}

# Print information for debug--------------------------
print_debug_info <- function(test, n1, n2, eval_results, low, high) {
    if (test == "intersection-union") {
        message(sprintf("Cluster size: %d, Clusters: %d, PMP1: %.3f, PMP2: %.3f, PMP3: %.3f, PMP4: %.3f, low: %d, high: %d",
                        n1, n2, eval_results$prop_PMP1, eval_results$prop_PMP2, 
                        eval_results$prop_PMP3, eval_results$prop_PMP4, low, high))
    } else if (test == "homogeneity") {
        message(sprintf("Cluster size: %d, Clusters: %d, BF12: %.3f, BF21: %.3f, low: %d, high: %d",
                        n1, n2, eval_results$prop_BF12, eval_results$prop_BF21, low, high))
    }
}

# Binary search: Updating sample size------------------------------------------
binary_search <- function(condition_met, test, fixed, n1, n2, low, high, max, eta,
                          current_eta, previous_eta, previous_high, min_sample){
    if (!condition_met) {
        message("Increasing sample size")
        if (fixed == "n1") {
            # Increase the number of clusters since eta is too small
            low <- n2                         #lower bound
            high <- high                      #higher bound
            n2 <- round2((low + high) / 2)     #point in the middle
            if (n2 %% 2 == 0) n2 <- n2 + 1 # To ensure number of clusters is even
            
            # Adjust higher bound when there is a ceiling effect
            if (low + n2 == high * 2) {
                low <- n2                         #lower bound
                if (previous_high > 0) {
                    high <- previous_high
                } else {
                    high <- max                       #higher bound
                }
                n2 <- round2((low + high) / 2)     #point in the middle
            }
            return(list(low = low,
                        high = high,
                        n1 = n1,
                        n2 = n2,
                        previous_eta = current_eta))
            
        } else if (fixed == "n2") {
                # Increase the cluster sizes since eta is too small
                low <- n1                        #lower bound
                high <- high                     #higher bound
                n1 <- round2((low + high) / 2)    #point in the middle
                
                # Adjust higher bound when there is a ceiling effect or no increase of power
                if ((low + n1 == high * 2) | (current_eta == previous_eta)) {
                    low <- n1                        #lower bound
                    #Set the higher bound based on the previous high or the maximum
                    if (previous_high > 0 ) {
                        high <- previous_high
                    } else {
                        high <- max
                    }
                    n1 <- round2((low + high) / 2)    # point in the middle
                }
                return(list(low = low,
                            high = high,
                            n1 = n1,
                            n2 = n2,
                            previous_eta = current_eta))
        }
    } else if (condition_met == TRUE) {
        previous_high <- high
        previous_eta <- current_eta
        return(list(previous_high = previous_high,
                    previous_eta = previous_eta))
    } # Finish condition met
}


# Final binary search---------------------------------------------------------
final_binary_search <- function(condition_met, test, fixed, n1, n2, low, high, max, eta,
                                current_eta, previous_eta, previous_high, min_sample){
    if (!condition_met) {
        message("Increasing sample size")
        if (fixed == "n1") {
            if (n2 == max)    { # If the sample size reaches the maximum
                return(list(
                    n1 = n1,
                    n2 = n2,
                    ultimate_sample_sizes = TRUE))
            } else {
                # Increase the number of clusters since eta is too small
                low <- n2                         #lower bound
                high <- high                      #higher bound
                n2 <- round2((low + high) / 2)     #point in the middle
                if (n2 %% 2 == 0) n2 <- n2 + 1 # To ensure number of clusters is even
                
                # Adjust higher bound when there is a ceiling effect
                if (low + n2 == high * 2) {
                    low <- n2                         #lower bound
                    if (previous_high > 0) {
                        high <- previous_high
                    } else {
                        high <- max                       #higher bound
                    }
                    n2 <- round2((low + high) / 2)     #point in the middle
                }
                return(list(low = low,
                            high = high,
                            n1 = n1,
                            n2 = n2,
                            ultimate_sample_sizes = FALSE,
                            previous_eta = current_eta))
            }
        } else if (fixed == "n2") {
            if (n1 == max)    {# If the sample size reaches the maximum
                return(list(low = min_sample,
                            high = max,
                            n1 = n1,
                            n2 = n2,
                            ultimate_sample_sizes = TRUE))
            } else {
                # Increase the cluster sizes since eta is too small
                low <- n1                        #lower bound
                high <- high                     #higher bound
                n1 <- round2((low + high) / 2)    #point in the middle
                
                # Adjust higher bound when there is a ceiling effect or no increase of power
                if ((low + n1 == high * 2) | (current_eta == previous_eta)) {
                    low <- n1                        #lower bound
                    #Set the higher bound based on the previous high or the maximum
                    if (previous_high > 0 ) {
                        high <- previous_high
                    } else {
                        high <- max
                    }
                    n1 <- round2((low + high) / 2)    # point in the middle
                }
                return(list(low = low,
                            high = high,
                            n1 = n1,
                            n2 = n2,
                            ultimate_sample_sizes = FALSE,
                            previous_eta = current_eta))
            }
        }
    } else if (condition_met == TRUE) {
        previous_high <- high
        print(c("previous:", previous_eta))
        
        if (fixed == "n1") {
            # Eta is close enough to the desired eta or
            # there is no change in eta and the lower bound is close to the middle point
            if ((current_eta - eta < 0.1 && n2 - low == 2) ||
                (previous_eta == current_eta && n2 - low == 2)){
                return(list(
                    n1 = n1,
                    n2 = n2,
                    ultimate_sample_sizes = TRUE))
                
            } else {
                # Decreasing to find the ultimate number of clusters
                message("Lowerign sample size")
                low <- low                         #lower bound
                high <- n2                         #higher bound
                n2 <- round2((low + high) / 2)      #point in the middle
                if (n2 %% 2 == 0) n2 <- n2 + 1
                if (n2 < 30) warning("The number of groups is less than 30.
                                             This may cause problems in convergence and singularity.")
                
                return(list(
                    low = low,
                    high = high,
                    n1 = n1,
                    n2 = n2,
                    ultimate_sample_sizes = FALSE,
                    previous_high = high, #The higher bound when the power criterion was met
                    previous_eta = current_eta
                ))
            }
        } else if (fixed == "n2") {
            # Eta is close enough to the desired eta or
            # there is no change in eta and the lower bound is close to the middle point or
            # reached the minimum number that meets the Bayesian power condition
            if ((current_eta - eta < 0.1 && n1 - low == 1) ||
                (current_eta == previous_eta && n1 - low == 1) ||
                (current_eta == previous_eta && low + n1 == high * 2)){
                return(list(low = min_sample,
                            high = max,
                            n1 = n1,
                            n2 = n2,
                            ultimate_sample_sizes = TRUE))
                
            } else {
                message("Lowerign sample size")
                # Decreasing the cluster size to find the ultimate sample size
                low <- low                         #lower bound
                high <- n1                         #higher bound
                n1 <- round2((low + high) / 2)      #point in the middle
                return(list(low = low,
                            high = high,
                            n1 = n1,
                            n2 = n2,
                            ultimate_sample_sizes = FALSE,
                            previous_high = high,
                            previous_eta = current_eta))
            }
        }
    } # Finish condition met
}
