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
    title <- "Final sample size"
    cat(paste("\n", title, "\n", sep = ""))
    row <- paste(rep("=", nchar(title)), collapse = "")
    cat(row, "\n")

    if (test == "intersection-union") {    # Print intersection-union
        cat("Sample Size Determination for Intersection-Union Test")
        print_hypotheses(list_hypo)
        results_df <- data.frame("n2" = integer(),
                                 "n1" = integer(),
                                 "P_BF12" = numeric(),
                                 "P_BF13" = numeric(),
                                 "P_BF14" = numeric(),
                                 "P_BF21" = numeric(),
                                 "P_BF31" = numeric(),
                                 "P_BF41" = numeric())

        names(results_df) <- c("Number of clusters", "Cluster size", 
                               paste("P(BF.12 >", object_result$BF_thres,"| H1)"),
                               paste("P(BF.13 >", object_result$BF_thres,"| H1)"),
                               paste("P(BF.14 >", object_result$BF_thres,"| H1)"),
                               paste("P(BF.21 >", object_result$BF_thres,"| H2)"),
                               paste("P(BF.31 >", object_result$BF_thres,"| H3)"),
                               paste("P(BF.41 >", object_result$BF_thres,"| H4)"))
        results_df[1, 1] <- object_result$n2 #n2
        results_df[1, 2] <- object_result$n1 #n1
        results_df[1, 3] <- object_result$Proportion.BF12 #BF_12
        results_df[1, 4] <- object_result$Proportion.BF13 #BF_13
        results_df[1, 5] <- object_result$Proportion.BF14 #BF_14
        results_df[1, 6] <- object_result$Proportion.BF21 #BF_21
        results_df[1, 7] <- object_result$Proportion.BF31 #BF_31
        results_df[1, 8] <- object_result$Proportion.BF41 #BF_41
        
        # Function to create a nicely formatted table with | as separators
        print_table <- function(df) {
            # Format the header row
            header <- paste(colnames(df), collapse = " | ")
            
            # Print the header row
            cat(header, "\n")
            
            # Print the separator line
            cat(strrep("-", nchar(header)), "\n")
            
            # Print each row with columns separated by |
            apply(df, 1, function(row) {
                cat(paste(row, collapse = " | "), "\n")
            })
        }
        
        # Function to split and print wide tables in chunks for better readability
        split_and_print_table <- function(df, chunk_size = 4) {
            # Split the data frame into smaller chunks by columns
            num_chunks <- ceiling(ncol(df) / chunk_size)
            
            for (i in 1:num_chunks) {
                # Get the current chunk of columns
                chunk <- df[, ((i-1)*chunk_size + 1):min(i*chunk_size, ncol(df)), drop = FALSE]
                print_table(chunk)
                cat("\n")  # Add a newline between chunks
            }
        }
        
        # Print the table (split into chunks if too wide)
        split_and_print_table(results_df, chunk_size = 4)
        
        # cat(knitr::kable(results_df, format = "pipe"))
        print(results_df)
        cat(paste(c("\u03B7 = ", object_result$eta)), "\n")
        cat("***********************************************************************", "\n")
    } else if (test == "homogeneity") {    # Print homogeneity of effect size
        cat("Hypotheses:", "\n")
        cat("    H1:", hypothesis1, "\n")
        cat("    Hc:", "complement", "\n")
        cat("Using a cluster size = ", object_result$n1, " and number of clusters = ", object_result$n2, "\n")
        cat("P (BF.1c > ", object_result$BF_thres, " | H1) = ", object_result$Proportion.BF1c, "\n")
        cat("On average the Posterior Model Probability is", mean(object_result$data_H1[, "PMP.1c"]))
        
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
