############################## PRINT RESULTS #################################

print_results <- function(object_result) {
    title <- "Final sample size"
    cat(paste("\n", title, "\n", sep = ""))
    row <- paste(rep("=", nchar(title)), collapse = "")
    cat(row, "\n")
    if (is.character(object_result[[length(object_result)]])) {   # Print for informative hypotheses
        cat("Hypotheses:", "\n")
        cat("    H1:", object_result[[5]][[1]], "\n")
        cat("    H2:", object_result[[5]][[2]], "\n")
        cat("Using cluster size = ", object_result$n1, " and number of clusters = ", object_result$n2, "\n")
        cat("P (BF.12 > ", object_result[[6]], " | H1) = ", object_result$Proportion.BF12, "\n")
    } else {                                    # Print for null vs informative
        n_object <- length(object_result)
        b_number <- n_object - 3
        results_matrix <- matrix(NA, nrow = b_number, ncol = 5)
        results_matrix[, 1] <- seq(b_number)
        object_result_b <- object_result[1:b_number]
        results_matrix[, 2] <- round(as.numeric(head(unlist(lapply(object_result_b, `[[`, 1)), b_number))) #n1
        results_matrix[, 3] <- round(as.numeric(head(unlist(lapply(object_result_b, `[[`, 2)), b_number))) #n2
        results_matrix[, 4] <- round(as.numeric(head(unlist(lapply(object_result_b, `[[`, 3)), b_number)), 3) #BF_01
        results_matrix[, 5] <- round(as.numeric(head(unlist(lapply(object_result_b, `[[`, 4)), b_number)), 3) #BF_10
        colnames(results_matrix) <- c("b", "n1", "n2", paste("P(BF.01 >", object_result[[(n_object - 1)]][[1]], "| H0) > ", object_result[[n_object]][[1]], sep = " "), 
                                      paste("P(BF.10 >", object_result[[(n_object - 1)]][[2]], "| H1) > ", object_result[[n_object]][[2]], sep = " "))
        
        cat("Hypotheses:", "\n")
        cat("    H0:", object_result[[b_number + 1]][[1]], "\n")
        cat("    H1:", object_result[[b_number + 1]][[2]], "\n")
        
        cat("***********************************************************************", "\n")
        print(format(results_matrix, justify = "centre"))
        cat("***********************************************************************", "\n")
        cat("n1: Cluster sizes", "\n")
        cat("n2: Number of clusters", "\n")
    }
}

# Test
# print_results(ssd_results_null)

# Print SSD multivariate CRT---------------------------------------------------
print_results_multiv <- function(object_result, test, list_hypo) {
    # title <- "Final sample size"
    # cat(paste("\n", title, "\n", sep = ""))
    # row <- paste(rep("=", nchar(title)), collapse = "")
    # cat(row, "\n")
    
    if (test == "intersection-union") {    # Print intersection-union
        title <- "Sample Size Determination for Intersection-Union Test"
        header <- paste(rep("*", nchar(title)), collapse = "")
        cat(header, "\n")
        cat(title,"\n")
        row <- paste(rep("-", nchar(title)), collapse = "")
        cat(row, "\n")
        print_hypotheses(list_hypo)
        
        # Create dataframe with results
        results_PMP <- data.frame("Hypothesis" = character(),
                                  "Prop_PMP" = numeric(),
                                  "error" = numeric())
        results_BF <- data.frame("h" = character(),
                                 "P_BF1m" = numeric(),
                                 "P_BFm1" = numeric())
        names(results_PMP) <- c("Hypothesis",
                                paste("P(PMP > ", object_result$pmp_thresh, "| H)"),
                                " Error")
        names(results_BF) <- c("Hypothesis", 
                               paste("BF.1m| H1"),
                               paste("BF.m1| Hm"))
        
        ## First column
        for (h in 2:length(list_hypo)) {
            results_BF[(h - 1), 1] <- paste0("m = ", h)
        }
        for (h in 1:length(list_hypo)) {
            results_PMP[h, 1] <- paste("Hypothesis", h)
        }
        
        ## Second column
        results_BF[1, 2] <- median(object_result$data$results_H1[, 2]) #BF_12
        results_BF[2, 2] <- median(object_result$data$results_H1[, 3]) #BF_13
        results_BF[3, 2] <- median(object_result$data$results_H1[, 4]) #BF_14
        
        results_PMP[1, 2] <- object_result$Proportion.PMP1
        results_PMP[2, 2] <- object_result$Proportion.PMP2
        results_PMP[3, 2] <- object_result$Proportion.PMP3
        results_PMP[4, 2] <- object_result$Proportion.PMP4
        
        ## Third column
        results_BF[1, 3] <- median(object_result$data$results_H2[, 2]) #BF_21
        results_BF[2, 3] <- median(object_result$data$results_H3[, 2]) #BF_31
        results_BF[3, 3] <- median(object_result$data$results_H4[, 2]) #BF_41
        
        results_PMP[1, 3] <- 1 - object_result$Proportion.PMP1
        results_PMP[2, 3] <- 1 - object_result$Proportion.PMP2
        results_PMP[3, 3] <- 1 - object_result$Proportion.PMP3
        results_PMP[4, 3] <- 1 - object_result$Proportion.PMP4
        
        # Print results
        cat("Number of clusters: ", object_result$n2, "\n")
        cat("Cluster size: ", object_result$n1, "\n")
        stargazer::stargazer(results_PMP, type = "text", summary = FALSE)
        # Print eta
        cat("\n", paste(c("\u03B7 = ", object_result$eta)), "\n")
        cat("\n")
        stargazer::stargazer(results_BF, type = "text", summary = FALSE)
        cat("***********************************************************************", "\n")
        # TODO: Add a conclusion of how to interpret the PMP and 1-PMP.
        
    } else if (test == "homogeneity") {    # Print homogeneity of effect size
        title <- "Sample Size Determination for Homogeneity of Effect Sizes Test"
        header <- paste(rep("*", nchar(title)), collapse = "")
        cat(header, "\n")
        cat(title, "\n")
        row <- paste(rep("-", nchar(title)), collapse = "")
        cat(row, "\n")
        print_hypotheses(list_hypo)
        
        # Create dataframe with results
        results_PMP <- data.frame("h" = character(),
                                  "Prop_PMP" = numeric(),
                                  "Error" = numeric())
        names(results_PMP) <- c("Hypothesis",
                                paste("P(PMP > )", object_result$pmp_thresh, "| H"),
                                "Error")
        results_BF <- data.frame("h" = character(),
                                 "Med_BF" = numeric())
        names(results_BF) <- c("Hypothesis",
                               "Median BF|H")
        
        # First columns
        for (h in 1:length(list_hypo)) {
            results_PMP[h, 1] <- paste("Hypothesis", h)
        }
        for (h in 1:length(list_hypo)) {
            results_BF[h, 1] <- paste("Hypothesis", h)
        }
        # Second columns
        results_PMP[1, 2] <- object_result$Proportion.PMP1
        results_PMP[2, 2] <- object_result$Proportion.PMP2
        
        results_BF[1, 2] <- median(object_result$data$results_H1[, 2]) #BF_12
        results_BF[2, 2] <- median(object_result$data$results_H2[, 2]) #BF_21
        
        # Third column
        results_PMP[1, 3] <- 1 - object_result$Proportion.PMP1
        results_PMP[2, 3] <- 1 - object_result$Proportion.PMP2
        
        # Print results
        cat("Number of clusters: ", object_result$n2, "\n")
        cat("Cluster size: ", object_result$n1, "\n")
        stargazer::stargazer(results_PMP, type = "text", summary = FALSE)
        # Print eta
        cat("\n", paste(c("\u03B7 = ", object_result$eta)), "\n")
        cat("\n")
        # Print BF
        stargazer::stargazer(results_BF, type = "text", summary = FALSE)
        cat("***********************************************************************", "\n")
        # TODO: Add a conclusion of how to interpret the PMP and 1-PMP.
        
    } else if (test == "omnibus") {                                    # Print for null vs informative
        title <- "Sample Size Determination for Omnibus Test"
        header <- paste(rep("*", nchar(title)), collapse = "")
        cat(header, "\n")
        cat(title,"\n")
        row <- paste(rep("-", nchar(title)), collapse = "")
        cat(row, "\n")
        print_hypotheses(list_hypo)
        
        # Create dataframe with results
        results_PMP <- data.frame("Hypothesis" = character(),
                                  "Prop_PMP" = numeric(),
                                  "error" = numeric())
        results_BF <- data.frame("h" = character(),
                                 "P_BF1m" = numeric(),
                                 "P_BFm1" = numeric())
        names(results_PMP) <- c("Hypothesis",
                                paste("P(PMP > )", object_result$pmp_thresh, "| H"),
                                " Error")
        names(results_BF) <- c("Hypothesis", 
                               paste("BF.m4| Hm"),
                               paste("BF.4m| H4"))
        
        ## First column
        for (h in 1:length(list_hypo)) {
            results_PMP[h, 1] <- paste("Hypothesis", h)
        }
        for (h in 1:(length(list_hypo) - 1)) {
            results_BF[(h - 1), 1] <- paste0("m = ", h)
        }
        
        ## Second column
        results_PMP[1, 2] <- object_result$Proportion.PMP1
        results_PMP[2, 2] <- object_result$Proportion.PMP2
        results_PMP[3, 2] <- object_result$Proportion.PMP3
        results_PMP[4, 2] <- object_result$Proportion.PMP4
        
        results_BF[1, 2] <- median(object_result$data$results_H1[, 2]) #BF_14
        results_BF[2, 2] <- median(object_result$data$results_H2[, 2]) #BF_24
        results_BF[3, 2] <- median(object_result$data$results_H3[, 2]) #BF_34
        
        ## Third column
        results_PMP[1, 3] <- 1 - object_result$Proportion.PMP1
        results_PMP[2, 3] <- 1 - object_result$Proportion.PMP2
        results_PMP[3, 3] <- 1 - object_result$Proportion.PMP3
        results_PMP[4, 3] <- 1 - object_result$Proportion.PMP4
        
        results_BF[1, 3] <- median(object_result$data$results_H4[, 2]) #BF_41
        results_BF[2, 3] <- median(object_result$data$results_H4[, 3]) #BF_42
        results_BF[3, 3] <- median(object_result$data$results_H4[, 4]) #BF_43
        
        # Print results
        cat("Number of clusters: ", object_result$n2, "\n")
        cat("Cluster size: ", object_result$n1, "\n")
        stargazer::stargazer(results_PMP, type = "text", summary = FALSE)
        # Print eta
        cat("\n", paste(c("\u03B7 = ", object_result$eta)), "\n")
        cat("\n")
        stargazer::stargazer(results_BF, type = "text", summary = FALSE)
        cat("***********************************************************************", "\n")
        # TODO: Add a conclusion of how to interpret the PMP and 1-PMP.
    }
}
