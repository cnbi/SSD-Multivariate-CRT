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
