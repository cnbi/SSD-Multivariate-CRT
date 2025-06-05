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
        # Extract the elements for outcome 1
        estimates_y1 <- data$theta$zeta[1:2]
        names(estimates_y1) <- c("intercept", "slope")
        
        # Extract the elements for outcome 2
        estimates_y2 <- data$theta$zeta[3:4]
        names(estimates_y2) <- c("Intercept", "slope")
        
        fixed_eff <- c(estimates_y1, estimates_y2)
        names(fixed_eff) <- c("Intercept1", "Slope1", "Intercept2", "Slope2")
    } else if (n_outcomes == 3) {
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
extract_var_cov <- function(estimations, n_outcomes) {
    if (n_outcomes == 2) {
        indexes_var_cov <- seq(2, 2*n_outcomes, by = 2)
        var_cov_matrix <- matrix(NA, n_outcomes, n_outcomes)
        for (i in seq(length(indexes_var_cov))) {
            for (j in seq(length(indexes_var_cov))) {
                row <- indexes_var_cov[i]
                col <- indexes_var_cov[j]
                var_cov_matrix[i, j] <- estimations$var_cov[row, col]
            }
        }
    } else if (n_outcomes == 3) {
        indexes_var_cov <- seq(2, 2*n_outcomes, by = 2)
        var_cov_matrix <- matrix(NA, n_outcomes, n_outcomes)
        for (i in seq(length(indexes_var_cov))) {
            for (j in seq(length(indexes_var_cov))) {
                row <- indexes_var_cov[i]
                col <- indexes_var_cov[j]
                var_cov_matrix[i, j] <- estimations$var_cov[row, col]
            }
        }
    }
    return(var_cov_matrix)
}


# Calculate empirical ICCs ---------------------------------
#
calc_ICCs <- function(estimations, n_outcomes) {
    emp_rho0 <- vector(mode = "numeric", length = n_outcomes)  # Empirical Outcome-specific ICC
    emp_rho1 <- matrix(NA, n_outcomes, n_outcomes)  # Empirical intersubject between-outcome ICC
    emp_rho2 <- matrix(NA, n_outcomes, n_outcomes)  # Empirical intrasubject between-outcome ICC
    
    # Calculate empirical ICCs
    emp_Sigma_phi <- estimations$theta$SigmaPhi
    emp_Sigma_e <- estimations$theta$SigmaE
    
    # Outcome-specific ICC
    for (y in 1:n_outcomes) {
        emp_rho0[y] <- emp_Sigma_phi[y ,y]/(emp_Sigma_phi[y ,y] + emp_Sigma_e[y, y])
    }
    if (n_outcomes == 2) {
        names(emp_rho0) <- c("y1", "y2")
    } else if (n_outcomes == 3) {
        names(emp_rho0) <- c("y1", "y2", "y3")
    }
    
    # Intersubject between-outcome ICC and Intrasubject between-outcome ICC
    for (y in 1:(n_outcomes - 1)) {
        for (y_prime in (y + 1):n_outcomes) {
            emp_rho1[y, y_prime] <- emp_Sigma_phi[y, y_prime]/(sqrt(emp_Sigma_phi[y, y] + emp_Sigma_e[y, y])*sqrt(emp_Sigma_phi[y_prime, y_prime] + emp_Sigma_e[y_prime, y_prime]))
            emp_rho1[y_prime, y] <- emp_rho1[y, y_prime]
            emp_rho2[y, y_prime] <- (emp_Sigma_phi[y, y_prime] + emp_Sigma_e[y, y_prime]) / (sqrt(emp_Sigma_phi[y, y] + emp_Sigma_e[y, y])*sqrt(emp_Sigma_phi[y_prime, y_prime] + emp_Sigma_e[y_prime, y_prime]))
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
        eigenvalues <- pmax(eigenvalues, 1e-6)  # Add a small constant (1e-6) to non-positive eigenvalues
        # Reconstruct the covariance matrix using the adjusted eigenvalues
        corrected_matrix <- eigenvectors %*% diag(eigenvalues) %*% t(eigenvectors)
        return(corrected_matrix)
    } else {
        return(cov_matrix)  # Return original matrix if it's already positive definite
    }
}