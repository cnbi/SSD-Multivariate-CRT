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

BF_multiv <- function(estimates, sigma, effective_n, hypothesis, pack, difference, test){
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
            Bf1u <- Bf$fit$Fit[1] / Bf$fit$Com[1]
            Bf2u <- Bf$fit$Fit[2] / Bf$fit$Com[2]
            Bf12 <- Bf1u / Bf2u
            Bf21 <- Bf2u / Bf1u
            Bf1c <- Bf$fit$BF.c[1]
            Bf2c <- Bf$fit$BF.c[2]
            PMP1c <- Bf$fit$PMPc[1]
            PMP2c <- Bf$fit$PMPc[2]
            
            if (test == "intersection-union") {
                Bf3u <- Bf$fit$Fit[3] / Bf$fit$Com[3]
                Bf4u <- Bf$fit$Fit[4] / Bf$fit$Com[4]
                Bf13 <- Bf1u / Bf3u
                Bf14 <- Bf1u / Bf4u
                Bf31 <- Bf3u / Bf1u
                Bf41 <- Bf4u / Bf1u
                Bf1c <- Bf$fit$BF.c[1]
                Bf3c <- Bf$fit$BF.c[3]
                Bf4c <- Bf$fit$BF.c[4]
                PMP3c <- Bf$fit$PMPc[3]
                PMP4c <- Bf$fit$PMPc[4]
            } else if (test == "omnibus") {
                Bf3u <- Bf$fit$Fit[3] / Bf$fit$Com[3]
                
            }
            
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
        if (test == "intersection-union") {
            
            if (any(is.na(c(Bf1u, Bf1c, Bf2u, Bf2c, Bf12, Bf21, Bf3u, Bf4u, Bf13, Bf14, Bf31, Bf41)))) {
                good_result <- FALSE
            } else if (any(is.nan(c(Bf1u, Bf1c, Bf2u, Bf2c, Bf12, Bf21, Bf3u, Bf4u, Bf13, Bf14, Bf31, Bf41)))) {
                good_result <- FALSE
            } else if (any(is.null(c(Bf1u, Bf1c, Bf2u, Bf2c, Bf12, Bf21, Bf3u, Bf4u, Bf13, Bf14, Bf31, Bf41)))) {
                good_result <- FALSE 
            } else {
                good_result <- TRUE
            }
        } else if (test == "homogeneity"){
            if (any(is.na(c(Bf1u, Bf1c, Bf2u, Bf2c, Bf12, Bf21)))) {
                good_result <- FALSE
            } else if (any(is.nan(c(Bf1u, Bf1c, Bf2u, Bf2c, Bf12, Bf21)))) {
                good_result <- FALSE
            } else if (any(is.null(c(Bf1u, Bf1c, Bf2u, Bf2c, Bf12, Bf21)))) {
                good_result <- FALSE 
            } else {
                good_result <- TRUE
            }
        }
    }

    if (test == "intersection-union") {
        results <- list(BF.1u = Bf1u, BF.2u = Bf2u, BF.3u = Bf3u, BF.4u = Bf4u, 
                        BF.12 = Bf12, BF.13 = Bf13, BF.14 = Bf14, 
                        BF.21 = Bf21, BF.31 = Bf31, BF.41 = Bf41,
                        BF.1c = Bf1c, BF.2c = Bf2c, BF.3c = Bf3c, BF.4c = Bf4c,
                        PMP.1c = PMP1c, PMP.2c = PMP2c, PMP.3c = PMP3c, PMP.4c = PMP4c)
    } else if (test == "homogeneity") {
        results <- list(BF.1u = Bf1u, BF.2u = Bf2u,
                        BF.12 = Bf12, BF.21 = Bf21,
                        BF.1c = Bf1c, BF.2c = Bf2c,
                        PMP.1c = PMP1c, PMP.2c = PMP2c)
    }

    return(results)
}
