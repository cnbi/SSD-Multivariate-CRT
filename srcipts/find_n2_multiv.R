############################ FIND SAMPLE SIZE ############################

#'effect_sizes: Vector with effect sizes (slopes).
#'n1: Common cluster sizes.
#'n2: Common number of clusters.
#'ndatasets: Number of data sets generated.
#'out_specific_ICC: Outcome/endpoint-specific intraclass correlation coefficient.
#'intersubj_between_outICC: Intersubject between-endpoint/outcome intraclass correlation coefficient.
#'intrasubj_between_outICC: Intrasubject between-endpoint/outcomes intraclass correlation coefficient.
#'Bayes_pack


SSD_mult_CRT <- function(test, effect_sizes, n1 = 15, n2 = 30, ndatasets = 1000, out_specific_ICC, 
                         intersubj_between_outICC, intrasubj_between_outICC,
                         BF_thresh = 3, eta = 0.8, fixed = "n1", difference = 0.2, max,
                         batch_size = 1000, seed, Bayes_pack) {
    # Libraries
    if (Bayes_pack == "bain") {
        if (!require("bain")) {install.packages("bain")}
    } else if (Bayes_pack == "BFpack") {
        if (!require("BFpack")) {install.packages("BFpack")}
    }
    if (!require("MASS")) {install.packages("MASS")} #Multivariate data generation
    if (!require("dplyr")) {install.packages("dplyr")} #?
    if (!require("numDeriv")) {install.packages("numDeriv")} #Hessian
    if (!require("nlme")) {install.packages("nlme")} #Multilevel model
    if (!require("lme4")) {install.packages("lme4")} #Multilevel model
    if (!require("stringr")) {install.packages("stringr")} #Multilevel model
    if (!requireNamespace("numDeriv", quietly = TRUE)) {
        install.packages("numDeriv")
    }
    if (!requireNamespace("nlme", quietly = TRUE)) {
        install.packages("nlme")
    }
    if (!requireNamespace("mvtnorm", quietly = TRUE)) {
        install.packages("mvtnorm")
    }
    
    # Functions
    source("srcipts/multiv_data_generation.R")
    source("srcipts/loglikelihood.R")
    source("srcipts/EM_algorithm.R")
    source("srcipts/aafbf.R")
    source("srcipts/helpers_functions.R")
    source("srcipts/print_results.R")
    
    
    # Warnings
    if (is.numeric(c(effect_sizes, n1, n2, ndatasets, out_specific_ICC, 
                     intersubj_between_outICC, intrasubj_between_outICC, BF_thresh, 
                     eta, max, batch_size)) == FALSE) 
        stop("All arguments, except 'fixed', must be numeric")
    # if (n2 %% 2 > 0) stop("The number of clusters must be even")
    #if (rho < 0) stop("The intraclass correlation must be a positive value")
    if (eta > 1) stop("The probability of exceeding Bayes Factor threshold cannot be larger than 1")
    if (eta < 0) stop("The probability of exceeding Bayes Factor threshold must be a positive value")
    if (is.character(fixed) == FALSE) stop("Fixed can only be a character indicating n1 or n2.")
    if (fixed %in% c("n1", "n2") == FALSE) stop("Fixed can only be a character indicating n1 or n2.")
    #if ((b == round(b)) == FALSE) stop("The fraction of information (b) must be an integer")
    #TODO: Add warnings for ICCs.
    
    # Hypotheses
    if (test == "intersection-union") {
        H1 <- "Slope1>0 & Slope2>0"
        H2 <- "Slope1>0 & Slope2<0"
        effect_sizesH2 <- effect_sizes * c(1, -1)
    } else if (test == "omnibus") {
        H1 <- "Slope1>0 & Slope2>0" #Change this
        H2 <- "Slope1>0 & Slope2<0" #Change this
    } else if (test == "homogeneity") {
        H1 <- hypothesis_maker(c("Slope1", "Slope2"), difference, "<")
        H2 <- hypothesis_maker(c("Slope1", "Slope2"), difference, ">")
    }
    
    #Starting values
    # Binary search start ------------------------------------------------------
    if (fixed == "n1") {
        min_sample <- 6                     # Minimum number of clusters
        low <- min_sample                   #lower bound
    } else if (fixed == "n2") {
        min_sample <- 5                     # Minimum cluster size
        low <- min_sample                   #lower bound
    }
    high <- max                    #higher bound
    
    # Staring values
    n_outcomes <- length(effect_sizes)
    # b <- 1
    previous_high <- 0
    previous_eta <- 0
    current_eta <- 0
    ultimate_sample_sizes <- FALSE
    
    results_H1 <- matrix(NA, ndatasets, 3)
    
    # Data generation
    while (ultimate_sample_sizes == FALSE) {
        # If H1 is true
        data_H1 <- do.call(gen_multiv_data, list(ndatasets, n1, n2, effect_sizes,
                                                 out_specific_ICC,
                                                 intersubj_between_outICC,
                                                 intrasubj_between_outICC,
                                                 n_outcomes, seed))
        
        # If H2 is true
        # data_H2 <- do.call(gen_multiv_data, list(ndatasets, n1, n2, effect_sizesH2, out_specific_ICC, intersubj_between_outICC, intrasubj_between_outICC,
        #                                          n_outcomes))
        
        # Effective sample size
        if (test == "intersection-union") {
            effective_n <- Map(effective_sample, list(n1), list(n2), data_H1$ICCs, list(n_outcomes))
            # effective_n <- effective_sample(n1, n2, data_H1$rho0, n_outcomes) #One for each outcome?
        } else if (test == "omnibus") {
            effective_n <- effective_sample(n1, n2)
        } else if (test == "homogeneity") {
            effective_n <- Map(effective_sample, list(n1), list(n2), data_H1$ICCs, list(n_outcomes), list(difference))
        }
        effective_n <- Map(min, effective_n)
        
        #Bayes factors------------------------------

        if (test == "intersection-union") {
            output_BF_H1 <- Map(BF_multiv, data_H1$estimations, data_H1$Sigma, effective_n, list(H1), list(Bayes_pack))
        } else if (test == "omnibus") {
            output_BF_H1 <- Map(BF_multiv, data_H1$estimations, data_H1$Sigma, effective_n, list(H1), list(Bayes_pack), list(difference))
        } else if ("homogeneity") {
            
            effective_n <- Map(effective_sample, list(n1), list(n2), data_H1$ICCs, list(n_outcomes), list(difference))
            }
        # output_BF_H1 <- BF_multiv(data_H1$estimations, data_H1$Sigma, effective_n, hypothesis = H1, pack = Bayes_pack)
        
        # Results ---------------------------------------------------------------------
        results_H1[, 1] <- unlist(lapply(output_BF_H1, extract_res, 1)) # Bayes factor H1vsHu
        results_H1[, 2] <- unlist(lapply(output_BF_H1, extract_res, 2)) # Bayes factor H1vsHc
        results_H1[, 3] <- unlist(lapply(output_BF_H1, extract_res, 3)) #posterior model probabilities of H1
        colnames(results_H1) <- c("BF.1u", "BF.1c", "PMP.1c")
        print("Bayes factor done!")
        # results_H0[, 1] <- unlist(lapply(output_AAFBF_H0, extract_res, 1)) # Bayes factor H1vsH0
        # results_H0[, 2] <- unlist(lapply(output_AAFBF_H0, extract_res, 4)) #posterior model probabilities of H1
        # results_H0[, 3] <- unlist(lapply(output_AAFBF_H0, extract_res, 2)) # Bayes factor H0vsH1
        # results_H0[, 4] <- unlist(lapply(output_AAFBF_H0, extract_res, 3)) #posterior model probabilities of H0
        # 
        # colnames(results_H0) <- c("BF.10", "PMP.1", "BF.01", "PMP.0")
        
        #Evaluation of condition -------------------------------------------
        # Proportion
        prop_BF1c <- length(which(results_H1[, "BF.1c"] > BF_thresh)) / ndatasets
        
        # Evaluation
        condition_met <- ifelse(prop_BF1c > eta, TRUE, FALSE)
        previous_eta <- current_eta
        current_eta <- prop_BF1c
        print("Bayes factor check!")

        # Binary search algorithm ------------------------------------------
        if (condition_met == FALSE) {
            print(c("Using cluster size:", n1, "and number of clusters:", n2,
                    "prop_BF1c: ", prop_BF1c, "low:", low, "high:", high))
            print("Increasing sample")
            if (fixed == "n1") {
                if (n2 == max)    { # If the sample size reaches the maximum
                    low <- min_sample
                    high <- max
                    ultimate_sample_sizes <- TRUE
                } else {
                    # Increase the number of clusters since eta is too small
                    low <- n2                         #lower bound
                    high <- high                      #higher bound
                    n2 <- round2((low + high) / 2)     #point in the middle
                    ifelse(n2 %% 2 == 0, n2 <- n2, n2 <- n2 + 1) # To ensure number of clusters is even
                    
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
                }
            } else if (fixed == "n2") {
                if (n1 == max)    {# If the sample size reaches the maximum
                    low <- min_sample
                    high <- max
                    ultimate_sample_sizes <- TRUE
                } else {
                    # Increase the cluster sizes since eta is too small
                    low <- n1                        #lower bound
                    high <- high                     #higher bound
                    n1 <- round2((low + high) / 2)    #point in the middle
                    
                    # Adjust higher bound when there is a ceiling effect
                    if ((low + n1 == high * 2) | (current_eta == previous_eta)) {
                        low <- n1                        #lower bound
                        #Set the higher bound based on the previous high or the maximum
                        if (previous_high > 0 ) {
                            high <- previous_high
                        } else {
                            high <- max
                        }
                        n1 <- round2((low + high) / 2)    #point in the middle
                    }
                }
            }
        } else if (condition_met == TRUE) {
            print(c("Using cluster size:", n1,
                    "and number of clusters:", n2,
                    "prop_BF1c: ", prop_BF1c,
                    "low: ", low, "high: ", high))
            previous_high <- high
            print("Lowerign sample")
            print(c("previous:", previous_eta))
            previous_eta <- current_eta
            
            if (fixed == "n1") {
                # Eta is close enough to the desired eta
                if (current_eta - eta < 0.1 && n2 - low == 2) {
                    low <- min_sample
                    high <- max
                    ultimate_sample_sizes <- TRUE
                    
                } else if (previous_eta == current_eta && n2 - low == 2) {
                    # If there is no change in eta and the lower bound is close to the middle point
                    low <- min_sample
                    high <- max
                    ultimate_sample_sizes <- TRUE
                    
                } else {
                    # Decreasing to find the ultimate number of clusters
                    low <- low                         #lower bound
                    high <- n2                         #higher bound
                    n2 <- round2((low + high) / 2)      #point in the middle
                    ifelse(n2 %% 2 == 0, n2 <- n2, n2 <- n2 + 1)
                    if (n2 < 30) warning("The number of groups is less than 30.
                                             This may cause problems in convergence and singularity.")
                    print("Lowering") # Eliminate later
                    
                }
            } else if (fixed == "n2") {
                # Eta is close enough to the desired eta
                if (current_eta - eta < 0.1 && n1 - low == 1) {
                    low <- min_sample
                    high <- max
                    ultimate_sample_sizes <- TRUE
                    
                } else if (current_eta == previous_eta && n1 - low == 1) {
                    # If there is no change in eta and the lower bound is close to the middle point
                    low <- min_sample
                    high <- max
                    ultimate_sample_sizes <- TRUE
                    
                } else if (current_eta == previous_eta && low + n2 == high * 2) {
                    # Reached the minimum number that meets the Bayesian power condition
                    low <- min_sample
                    high <- max
                    ultimate_sample_sizes <- TRUE
                    
                } else {
                    # Decreasing the cluster size to find the ultimate sample size
                    low <- low                         #lower bound
                    high <- n1                         #higher bound
                    n1 <- round2((low + high) / 2)      #point in the middle
                    print("Lowering") # Eliminate later
                }
            }
        } # Finish condition met
        
    } # Finish while loop ultimate_sample_size

    SSD_object <- list("n1" = n1,
                       "n2" = n2,
                       "Proportion.BF1c" = prop_BF1c,
                       "BF_thres" = BF_thresh,
                       "eta" = eta,
                       "data_H1" = results_H1)
    
    # Final output
    print_results_multiv(SSD_object, test, H1)
    invisible(SSD_object)
}
