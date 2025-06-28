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
                         seed, Bayes_pack) {
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
    source("scripts/multiv_data_generation.R")
    source("scripts/loglikelihood.R")
    source("scripts/EM_algorithm.R")
    source("scripts/aafbf.R")
    source("scripts/helpers_functions.R")
    source("scripts/print_results.R")
    
    
    # Warnings
    if (is.numeric(c(effect_sizes, n1, n2, ndatasets, out_specific_ICC, 
                     intersubj_between_outICC, intrasubj_between_outICC, BF_thresh, 
                     eta, max)) == FALSE) 
        stop("All arguments, except 'fixed', must be numeric")
    # if (n2 %% 2 > 0) stop("The number of clusters must be even")
    #if (rho < 0) stop("The intraclass correlation must be a positive value")
    if (eta > 1) stop("The probability of exceeding Bayes Factor threshold cannot be larger than 1")
    if (eta < 0) stop("The probability of exceeding Bayes Factor threshold must be a positive value")
    if (is.character(fixed) == FALSE) stop("Fixed can only be a character indicating n1 or n2.")
    if (fixed %in% c("n1", "n2") == FALSE) stop("Fixed can only be a character indicating n1 or n2.")
    if ((test == "homogeneity") && (effect_sizes[1] < effect_sizes[2]))
        stop("Effect size 1 must be larger than effect size 2")
    
    #TODO: Add warnings for ICCs.
    
    # Hypotheses
    if (test == "intersection-union") {
        H1 <- "Slope1>0 & Slope2>0"
        H2 <- "Slope1>0 & Slope2<0"
        H3 <- "Slope1<0 & Slope2>0"
        H4 <- "Slope1<0 & Slope2<0"
        effect_sizesH2 <- effect_sizes * c(1, -1)
        effect_sizesH3 <- effect_sizes * c(-1, 1)
        effect_sizesH4 <- effect_sizes * c(-1, -1)
    } else if (test == "omnibus") {
        H1 <- "Slope1>0 & Slope2>0" #Change this
        H2 <- "Slope1>0 & Slope2<0" #Change this
    } else if (test == "homogeneity") {
        H1 <- hypothesis_maker(c("Slope1", "Slope2"), difference, "<")
        H2 <- hypothesis_maker(c("Slope1", "Slope2"), difference, ">")
        effect_sizesH2 <- c(effect_sizes[1] + difference, effect_sizes[2])
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
    
    if (test == "intersection-union") {
        results_H1 <- matrix(NA, ndatasets, 5)
        results_H2 <- matrix(NA, ndatasets, 3)
        results_H3 <- matrix(NA, ndatasets, 3)
        results_H4 <- matrix(NA, ndatasets, 3)
    } else if (test == "homogeneity") {
        results_H1 <- matrix(NA, ndatasets, 3)
        results_H2 <- matrix(NA, ndatasets, 3)
    }
    
    seed <- ifelse(missing(seed), 1920, seed)
    
    
    # Data generation
    while (ultimate_sample_sizes == FALSE) {
        if (test == "intersection-union") {
            # If H1 is true
            data_H1 <- do.call(gen_multiv_data, list(ndatasets, n1, n2, effect_sizes,
                                                     out_specific_ICC,
                                                     intersubj_between_outICC,
                                                     intrasubj_between_outICC,
                                                     n_outcomes, seed))
            
            # If H2 is true
            data_H2 <- do.call(gen_multiv_data, list(ndatasets, n1, n2, effect_sizesH2, 
                                                     out_specific_ICC, 
                                                     intersubj_between_outICC, 
                                                     intrasubj_between_outICC,
                                                     n_outcomes, seed))
            
            data_H3 <- do.call(gen_multiv_data, list(ndatasets, n1, n2, effect_sizesH3, 
                                                     out_specific_ICC, 
                                                     intersubj_between_outICC, 
                                                     intrasubj_between_outICC,
                                                     n_outcomes, seed))
            
            data_H4 <- do.call(gen_multiv_data, list(ndatasets, n1, n2, effect_sizesH4, 
                                                     out_specific_ICC, 
                                                     intersubj_between_outICC, 
                                                     intrasubj_between_outICC,
                                                     n_outcomes, seed))
            
        } else if (test == "homogeneity") {
            data_H1 <- do.call(gen_multiv_data, list(ndatasets, n1, n2, effect_sizes,
                                                     out_specific_ICC,
                                                     intersubj_between_outICC,
                                                     intrasubj_between_outICC,
                                                     n_outcomes, seed))
            
            # If H2 is true
            data_H2 <- do.call(gen_multiv_data, list(ndatasets, n1, n2, effect_sizesH2, 
                                                     out_specific_ICC, 
                                                     intersubj_between_outICC, 
                                                     intrasubj_between_outICC,
                                                     n_outcomes, seed))
            
        }
        
        # Effective sample size
        if (test == "intersection-union") {
            effective_nH1 <- Map(effective_sample, list(n1), list(n2), data_H1$ICCs, list(n_outcomes))
            effective_nH1 <- Map(min, effective_nH1)
            effective_nH2 <- Map(effective_sample, list(n1), list(n2), data_H2$ICCs, list(n_outcomes))
            effective_nH2 <- Map(min, effective_nH2)
            effective_nH3 <- Map(effective_sample, list(n1), list(n2), data_H3$ICCs, list(n_outcomes))
            effective_nH3 <- Map(min, effective_nH3)
            effective_nH4 <- Map(effective_sample, list(n1), list(n2), data_H4$ICCs, list(n_outcomes))
            effective_nH4 <- Map(min, effective_nH4)
            
        } else if (test == "omnibus") {
            effective_n <- effective_sample(n1, n2)
            effective_n <- Map(min, effective_n)
            
        } else if (test == "homogeneity") {
            effective_nH1 <- Map(effective_sample, list(n1), list(n2), data_H1$ICCs, list(n_outcomes))
            effective_nH1 <- Map(min, effective_nH1)
            effective_nH2 <- Map(effective_sample, list(n1), list(n2), data_H2$ICCs, list(n_outcomes))
            effective_nH2 <- Map(min, effective_nH2)
        }
        print("Effective sample size done")
        
        #Bayes factors-----------------------------------------------------------------
        
        if (test == "intersection-union") {
            output_BF_H1 <- Map(BF_multiv, data_H1$estimations, data_H1$Sigma, effective_nH1, list(paste0(H1, ";", H2, ";", H3, ";", H4)), list(Bayes_pack), list(test))
            output_BF_H2 <- Map(BF_multiv, data_H2$estimations, data_H2$Sigma, effective_nH2, list(paste0(H1, ";", H2, ";", H3, ";", H4)), list(Bayes_pack), list(test))
            output_BF_H3 <- Map(BF_multiv, data_H3$estimations, data_H3$Sigma, effective_nH3, list(paste0(H1, ";", H2, ";", H3, ";", H4)), list(Bayes_pack), list(test))
            output_BF_H4 <- Map(BF_multiv, data_H4$estimations, data_H4$Sigma, effective_nH4, list(paste0(H1, ";", H2, ";", H3, ";", H4)), list(Bayes_pack), list(test))
            
        } else if (test == "omnibus") {
            output_BF_H1 <- Map(BF_multiv, data_H1$estimations, data_H1$Sigma, effective_n, list(H1), list(Bayes_pack))
            
        } else if (test == "homogeneity") {
            H1_test <- paste0("Slope1-Slope2<", as.character(difference + 0.001)) #The hypothesis includes when the difference is exactly the difference specified
            H2_test <- paste0("Slope1-Slope2>", as.character(difference + 0.001)) #The hypothesis includes when the difference is exactly the difference specified
            output_BF_H1 <- Map(BF_multiv, data_H1$estimations, data_H1$Sigma, effective_nH1, list(paste0(H1_test, ";", H2_test)), list(Bayes_pack), list(difference), list(test))
            output_BF_H2 <- Map(BF_multiv, data_H2$estimations, data_H2$Sigma, effective_nH2, list(paste0(H1_test, ";", H2_test)), list(Bayes_pack), list(difference), list(test))
        }
        
        # Results ---------------------------------------------------------------------
        if (test == "intersection-union") {
            results_H1[, 1] <- unlist(lapply(output_BF_H1, extract_res, 1)) # Bayes factor H1vsHu
            results_H1[, 2] <- unlist(lapply(output_BF_H1, extract_res, 5)) # Bayes factor H1vsH2
            results_H1[, 3] <- unlist(lapply(output_BF_H1, extract_res, 6)) # Bayes factor H1vsH3
            results_H1[, 4] <- unlist(lapply(output_BF_H1, extract_res, 7)) # Bayes factor H1vsH4
            results_H1[, 5] <- unlist(lapply(output_BF_H1, extract_res, 11)) # Bayes factor H1vsHc
            #results_H1[, 6] <- unlist(lapply(output_BF_H1, extract_res, 15)) #posterior model probabilities of H1vsHc
            colnames(results_H1) <- c("BF.1u", "BF.12", "BF.13", "BF.14", "BF.1c")
            
            results_H2[, 1] <- unlist(lapply(output_BF_H2, extract_res, 2)) # Bayes factor H2vsHu
            results_H2[, 2] <- unlist(lapply(output_BF_H2, extract_res, 8)) # Bayes factor H2vsH1
            results_H2[, 3] <- unlist(lapply(output_BF_H2, extract_res, 12)) # Bayes factor H2vsHc
            #results_H2[, 4] <- unlist(lapply(output_BF_H2, extract_res, 16)) #posterior model probabilities of H2vsHc
            colnames(results_H2) <- c("BF.2u", "BF.21", "BF.2c")
            
            results_H3[, 1] <- unlist(lapply(output_BF_H3, extract_res, 3)) # Bayes factor H3vsHu
            results_H3[, 2] <- unlist(lapply(output_BF_H3, extract_res, 9)) # Bayes factor H3vsH1
            results_H3[, 3] <- unlist(lapply(output_BF_H3, extract_res, 13)) # Bayes factor H3vsHc
            #results_H3[, 4] <- unlist(lapply(output_BF_H3, extract_res, 17)) #posterior model probabilities of H3vsHc
            colnames(results_H3) <- c("BF.3u", "BF.31", "BF.3c")
            
            results_H4[, 1] <- unlist(lapply(output_BF_H4, extract_res, 4)) # Bayes factor H4vsHu
            results_H4[, 2] <- unlist(lapply(output_BF_H4, extract_res, 10)) # Bayes factor H4vsH1
            results_H4[, 3] <- unlist(lapply(output_BF_H4, extract_res, 14)) # Bayes factor H4vsHc
            #results_H4[, 4] <- unlist(lapply(output_BF_H4, extract_res, 18)) #posterior model probabilities of H4vsHc
            colnames(results_H4) <- c("BF.4u", "BF.41", "BF.4c")
            
        } else if (test == "homogeneity") {
            results_H1[, 1] <- unlist(lapply(output_BF_H1, extract_res, 1)) # Bayes factor H1vsHu
            results_H1[, 2] <- unlist(lapply(output_BF_H1, extract_res, 3)) # Bayes factor H1vsH2
            results_H1[, 3] <- unlist(lapply(output_BF_H1, extract_res, 5)) # Bayes factor H1vsHc
            colnames(results_H1) <- c("BF.1u", "BF.12", "BF.1c")
            
            results_H2[, 1] <- unlist(lapply(output_BF_H2, extract_res, 2)) # Bayes factor H2vsHu
            results_H2[, 2] <- unlist(lapply(output_BF_H2, extract_res, 4)) # Bayes factor H2vsH1
            results_H2[, 3] <- unlist(lapply(output_BF_H2, extract_res, 6)) # Bayes factor H2vsHc
            colnames(results_H2) <- c("BF.2u", "BF.21", "BF.2c")
            
        }
        
        print("Bayes factor done!")
        
        #Evaluation of condition -------------------------------------------
        # Proportion
        if (test == "intersection-union") {
            prop_BF12 <- length(which(results_H1[, "BF.12"] > BF_thresh)) / ndatasets
            prop_BF13 <- length(which(results_H1[, "BF.13"] > BF_thresh)) / ndatasets
            prop_BF14 <- length(which(results_H1[, "BF.14"] > BF_thresh)) / ndatasets
            prop_BF21 <- length(which(results_H2[, "BF.21"] > BF_thresh)) / ndatasets
            prop_BF31 <- length(which(results_H3[, "BF.31"] > BF_thresh)) / ndatasets
            prop_BF41 <- length(which(results_H4[, "BF.41"] > BF_thresh)) / ndatasets
        } else if (test == "homogeneity") {
            prop_BF12 <- length(which(results_H1[, "BF.12"] > BF_thresh)) / ndatasets
            prop_BF21 <- length(which(results_H2[, "BF.21"] > BF_thresh)) / ndatasets
        }
        
        # Evaluation
        if (test == "intersection-union") {
            condition_met <- ifelse(prop_BF12 > eta & prop_BF13 > eta & prop_BF14 > eta & prop_BF21 > eta & prop_BF31 > eta & prop_BF41 > eta, TRUE, FALSE)
            previous_eta <- current_eta
            if (prop_BF12 < eta & prop_BF21 < eta & prop_BF31 < eta & prop_BF41 < eta) {
                current_eta <- min(prop_BF12, prop_BF21, prop_BF31, prop_BF41, prop_BF13, prop_BF14)
            } else if (prop_BF12 < eta | prop_BF21 < eta | prop_BF31 < eta | prop_BF41 < eta | prop_BF13 | prop_BF14) {
                current_eta <- min(prop_BF12, prop_BF21, prop_BF31, prop_BF41, prop_BF13, prop_BF14)
            } else if (condition_met) {
                current_eta <- min(prop_BF12, prop_BF21, prop_BF31, prop_BF41, prop_BF13, prop_BF14)
            }
        } else if (test == "homogeneity") {
            #TODO: Change these
            condition_met <- ifelse(prop_BF12 > eta & prop_BF21 > eta, TRUE, FALSE)
            previous_eta <- current_eta
            if (prop_BF12 < eta & prop_BF21 < eta) {
                current_eta <- min(prop_BF12, prop_BF21)
            } else if (prop_BF12 < eta | prop_BF21 < eta) {
                current_eta <- min(prop_BF12, prop_BF21)
            } else if (condition_met) {
                current_eta <- min(prop_BF12, prop_BF21)
            }
        }
        print("Bayes factor check!")
        
        # Binary search algorithm ------------------------------------------
        if (condition_met == FALSE) {
            if (test == "intersection-union") {
                print(c("Using cluster size:", n1, 
                        "and number of clusters:", n2,
                        "prop_BF12: ", prop_BF12,
                        "prop_BF13: ", prop_BF13,
                        "prop_BF14: ", prop_BF14,
                        "prop_BF21: ", prop_BF21,
                        "prop_BF31: ", prop_BF31,
                        "prop_BF41: ", prop_BF41,
                        "low:", low, "high:", high))
            } else if (test == "homogeneity") {
                print(c("Using cluster size:", n1, 
                        "and number of clusters:", n2,
                        "prop_BF12: ", prop_BF12,
                        "prop_BF21: ", prop_BF21,
                        "low:", low, "high:", high))
            }
            
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
            if (test == "intersection-union") {
                print(c("Using cluster size:", n1,
                        "and number of clusters:", n2,
                        "prop_BF12: ", prop_BF12,
                        "prop_BF13: ", prop_BF13,
                        "prop_BF14: ", prop_BF14,
                        "prop_BF21: ", prop_BF21,
                        "prop_BF31: ", prop_BF31,
                        "prop_BF41: ", prop_BF41,
                        "low: ", low, "high: ", high))
            } else if (test == "homogeneity") {
                print(c("Using cluster size:", n1, 
                        "and number of clusters:", n2,
                        "prop_BF12: ", prop_BF12,
                        "prop_BF21: ", prop_BF21,
                        "low:", low, "high:", high))
            }
            
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
    
    if (test == "intersection-union") {
        SSD_object <- list("n1" = n1,
                           "n2" = n2,
                           "Proportion.BF12" = prop_BF12,
                           "Proportion.BF13" = prop_BF13,
                           "Proportion.BF14" = prop_BF14,
                           "Proportion.BF21" = prop_BF21,
                           "Proportion.BF31" = prop_BF31,
                           "Proportion.BF41" = prop_BF41,
                           "BF_thres" = BF_thresh,
                           "eta" = eta,
                           "data" = list(results_H1 = results_H1,
                                         results_H2 = results_H2,
                                         results_H3 = results_H3,
                                         results_H4 = results_H4))
        
    } else if (test == "homogeneity") {
        SSD_object <- list("n1" = n1,
                           "n2" = n2,
                           "Proportion.BF12" = prop_BF12,
                           "Proportion.BF21" = prop_BF21,
                           "BF_thres" = BF_thresh,
                           "eta" = eta,
                           "data" = list(results_H1 = results_H1,
                                         results_H2 = results_H2))
    }
    
    # Final output
    if (test == "intersection-union") {
        print_results_multiv(SSD_object, test, list(H1, H2, H3, H4))
    } else if (test == "homogeneity") {
        print_results_multiv(SSD_object, test, list(H1, H2))
    }
    invisible(SSD_object)
}
