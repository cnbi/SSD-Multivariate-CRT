####################### CALCULATE AAFBF #####################################

# type: String that indicates the type of comparison of hypotheses. "equality" 
#       test for only equality vs. only inequality, "inequalities" test for only
#       inequality vs. only inequality.
# estimates: R object with estimates.
# n: numeric. Effective sample size.
# sigma: list with covariances.
# b: Numerical. Fraction of information to specify the prior.
# n_eff: Effective sample size.

calc_aafbf <- function(type, estimates, sigma, b, n_eff, outcome_type) {
    if (type == "Inequalities") {
        # Complexities
        comp1 <- .5
        comp2 <- .5
        
        # Fit
        browser()
        fit2 <- pnorm(0, mean = estimates[1], sd = sqrt(sigma[[1]]))
        fit1 <- 1 - fit2
        
        # Calculation BFs
        AAFBF1u <- (fit1 / comp1)
        AAFBF2u <- (fit2 / comp2)
        AAFBF12 <- AAFBF1u / AAFBF2u
        AAFBF21 <- 1 / AAFBF12
        
        #Calculation of PMPs
        pmp1 <- AAFBF1u / (AAFBF1u + AAFBF2u)
        pmp2 <- 1 - pmp1
        output <- list(bf.12 = AAFBF12, bf.21 = AAFBF21, pmp1 = pmp1, pmp2 = pmp2)
        
    } else if (type == "Equality") {
        b_calc <- b * 1 / n_eff                   # Calculate b
        # Complexities
        comp0 <- dnorm(0, mean = 0, sd = sqrt(sigma[[1]] / b_calc))     # overlap of parameter under H0 and unconstrained prior -> density of the prior under Hu at the focal point 0
        comp1 <- 1 - pnorm(0, mean = 0, sd = sqrt(sigma[[2]] / b_calc))
        
        # Fit
        fit0 <- dnorm(0, mean = estimates[[1]], sd = sqrt(sigma[[1]])) # overlap of parameter under H0 and posterior -> density of the posterior at focal point 0
        fit1 <- 1 - pnorm(0, mean = estimates[[2]], sd = sqrt(sigma[[2]])) # the fit is equal to 1 - the fit of the complement
        
        # Calculation of BFs
        AAFBF0u <- fit0 / comp0                    # AAFBF of H0 vs Hu
        AAFBF1u <- fit1 / comp1                    # AAFBF of H1 vs Hu
        AAFBF01 <- AAFBF0u / AAFBF1u
        AAFBF10 <- 1 / AAFBF01
        
        # Calculation of PMPs
        pmp0 <- AAFBF0u / (AAFBF0u + AAFBF1u)
        pmp1 <- 1 - pmp0
        output <- list(bf.10 = AAFBF10, bf.01 = AAFBF01, pmp0 = pmp0, pmp1 = pmp1)
        
        
        
    }
    #Output
    return(output)
}


############## MULTIVARIATE TESTING ########################

# estimates: Vector numeric with the estimates
# sigma: Matrix with the variances and covariances

BF_multiv <- function(estimates, sigma, effective_n, hypothesis, pack, difference){
    name_parameters <- unique(stringr::str_extract_all(hypothesis, "\\b[[:alnum:]_]+\\b")[[1]])
    estimates <- estimates[names(estimates) %in% name_parameters]
    good_result <- FALSE
    while (good_result == FALSE) {
        if (pack == "BFpack") {
            # Using BFpack
            Bf <- BF(estimates, Sigma = sigma, n = effective_n, hypothesis = hypothesis)
            Bf1u <- Bf$BFtable_confirmatory[1, 6]
            Bf1c <- Bf$BFmatrix_confirmatory[1, 2]
            Bf_1c <- Bf$BFtable_confirmatory[1, 6] / Bf$BFtable_confirmatory[2, 6]
            PMP <- Bf$BFtable_confirmatory[1, 8]
        } else if (pack == "bain") {
            # Using bain
            Bf <- bain(estimates, hypothesis, n = effective_n, Sigma = sigma)
            Bf1u <- Bf$fit$Fit[1]/Bf$fit$Com[1]
            Bf1c <- Bf$fit$BF.c[1]
            Bf_1c <- (Bf$fit$Fit[1]/Bf$fit$Com[1])/(Bf$fit$Fit[3]/Bf$fit$Com[3])
            PMP <- Bf$fit$PMPc[1]
        } else {
            # Using my own code
            # # Complexities
            # complexity_h1 <- 1 - pmvnorm(lower = c(0, 0), upper = c(Inf, Inf), mean = c(0, 0), sd = sqrt(sigma[[2]] / b_calc))
            # 
            # # Fit
            # fit_h1 <- 1 - pmvnorm(lower = c(0, 0), upper = c(Inf, Inf), mean = estimates, sigma = sigma) # i could include lower as an argument that can be changed
            # fit_hc <- pmvnorm(lower = c(-Inf, -Inf), upper = c(0, 0), mean = estimates, sigma = sigma)
            # 
            # # Calculation of BFs
            # 
            # 
            # # Calculataion of PMPs
        }
        if (any(is.na(c(Bf1u, Bf1c, PMP)))) {
            good_result <- FALSE
        } else if (any(is.nan(c(Bf1u, Bf1c, PMP)))) {
            good_result <- FALSE
        } else if (any(is.null(c(Bf1u, Bf1c, PMP)))) {
            good_result <- FALSE
        } else {good_result <- TRUE}
    }
    results <- list(BF.1u = Bf1u, BF.1c = Bf1c, PMP.1c = PMP)
    return(results)
}
