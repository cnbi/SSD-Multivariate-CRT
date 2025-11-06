############# TESTS###########################################


source("find_n2_multiv.R")
options(error = recover)

# Intersection-union test ------------------------------------------------------
# Hypothesis true and bain
start_time <- Sys.time()
SSD_mult_CRT("intersection-union",
             effect_sizes = c(0.2, 0.5),
             n1 = 10,
             n2 = 26,
             ndatasets = 5,
             out_specific_ICC = 0.1,
             intersubj_between_outICC = 0.08,
             intrasubj_between_outICC = 0.05,
             BF_thresh = 3,
             eta = 0.8,
             fixed = "n1",
             max = 50,
             Bayes_pack = "bain")
end_time <- Sys.time()
as.numeric(difftime(end_time, start_time, units = "mins"))

# Low effect size

SSD_mult_CRT("intersection-union",
             effect_sizes = c(0.2, 0.2),
             n1 = 10,
             n2 = 26,
             ndatasets = 5,
             out_specific_ICC = 0.1,
             intersubj_between_outICC = 0.08,
             intrasubj_between_outICC = 0.25,
             BF_thresh = 3,
             eta = 0.8,
             fixed = "n1",
             max = 50,
             Bayes_pack = "bain")

# Is the order in the outcomes affect the sample size?
start_time <- Sys.time()
a <- SSD_mult_CRT("intersection-union",
             effect_sizes = c(0.2, 0.4),
             n1 = 10,
             n2 = 26,
             ndatasets = 5,
             out_specific_ICC = 0.1,
             intersubj_between_outICC = 0.08,
             intrasubj_between_outICC = 0.05,
             BF_thresh = 3,
             eta = 0.8,
             fixed = "n1",
             max = 50,
             Bayes_pack = "bain")
end_time <- Sys.time()
as.numeric(difftime(end_time, start_time, units = "mins"))

# Medium effect sizes
start_time <- Sys.time()
b <- SSD_mult_CRT("intersection-union",
             effect_sizes = c(0.4, 0.2),
             n1 = 10,
             n2 = 26,
             ndatasets = 5,
             out_specific_ICC = 0.1,
             intersubj_between_outICC = 0.08,
             intrasubj_between_outICC = 0.05,
             BF_thresh = 3,
             eta = 0.8,
             fixed = "n1",
             max = 50,
             Bayes_pack = "bain")
end_time <- Sys.time()
as.numeric(difftime(end_time, start_time, units = "mins"))


start_time <- Sys.time()
a <- SSD_mult_CRT("intersection-union",
                  effect_sizes = c(0.25, 0.7),
                  n1 = 10,
                  n2 = 26,
                  ndatasets = 5,
                  out_specific_ICC = 0.1,
                  intersubj_between_outICC = 0.08,
                  intrasubj_between_outICC = 0.05,
                  BF_thresh = 3,
                  eta = 0.8,
                  fixed = "n1",
                  max = 50,
                  Bayes_pack = "bain")
end_time <- Sys.time()
as.numeric(difftime(end_time, start_time, units = "mins"))


start_time <- Sys.time()
a <- SSD_mult_CRT("intersection-union",
                  effect_sizes = c(0.33, 0.66),
                  n1 = 10,
                  n2 = 26,
                  ndatasets = 5,
                  out_specific_ICC = 0.1,
                  intersubj_between_outICC = 0.08,
                  intrasubj_between_outICC = 0.05,
                  BF_thresh = 3,
                  eta = 0.8,
                  fixed = "n1",
                  max = 50,
                  Bayes_pack = "bain")
end_time <- Sys.time()
as.numeric(difftime(end_time, start_time, units = "mins"))

start_time <- Sys.time()
b <- SSD_mult_CRT("intersection-union",
                  effect_sizes = c(0.66, 0.33),
                  n1 = 10,
                  n2 = 26,
                  ndatasets = 5,
                  out_specific_ICC = 0.1,
                  intersubj_between_outICC = 0.08,
                  intrasubj_between_outICC = 0.05,
                  BF_thresh = 3,
                  eta = 0.8,
                  fixed = "n1",
                  max = 50,
                  Bayes_pack = "bain")
end_time <- Sys.time()
as.numeric(difftime(end_time, start_time, units = "mins"))

set.seed(2916)
test_random <- round(runif(n = 5, min = 1, max = nrow(design_matrix_n2)))
results <- vector("list", 5)
index <- 1
for (i in test_random) {
    results[[index]] <- SSD_mult_CRT("intersection-union",
                 effect_sizes = c(design_matrix_n2[i, 1], design_matrix_n2[i, 2]),
                 n1 = design_matrix_n2[i, 9],
                 n2 = 26,
                 ndatasets = 100,
                 out_specific_ICC = c(design_matrix_n2[i, 3], 0.1),
                 intersubj_between_outICC = design_matrix_n2[i, 4],
                 intrasubj_between_outICC = design_matrix_n2[i, 5],
                 pmp_thresh = design_matrix_n2[i, 6],
                 eta = design_matrix_n2[i, 7],
                 fixed = as.character(design_matrix_n2[i, 8]),
                 max = 100,
                 master.seed = 1629,
                 Bayes_pack = "bain")
    index <- index + 1
}

# Homogeneity of effect size test ----------------------------------------------
# Testing warning
a <- SSD_mult_CRT("homogeneity", effect_sizes = c(0.4, 0.5),
                  n1 = 10, n2 = 30, ndatasets = 10, out_specific_ICC = 0.025, 
                  intersubj_between_outICC = 0.022, intrasubj_between_outICC = 0.5,
                  BF_thresh = 5, eta = 0.8, fixed = "n1", difference = 0.2,
                  max = 100, batch_size = 100, seed = 1855, Bayes_pack = "bain")

# H1 must be true: the effect sizes are similar
a <- SSD_mult_CRT("homogeneity", effect_sizes = c(0.5, 0.4),
                  n1 = 10, n2 = 30, ndatasets = 10, out_specific_ICC = 0.025, 
                  intersubj_between_outICC = 0.022, intrasubj_between_outICC = 0.5,
                  BF_thresh = 5, eta = 0.8, fixed = "n1", difference = 0.2,
                  max = 500, batch_size = 100, seed = 1855, Bayes_pack = "bain")


# What happens if I have exactly the same difference as the one specified?
a <- SSD_mult_CRT("homogeneity", effect_sizes = c(0.5, 0.3),
                  n1 = 10, n2 = 30, ndatasets = 10, out_specific_ICC = 0.025, 
                  intersubj_between_outICC = 0.022, intrasubj_between_outICC = 0.5,
                  BF_thresh = 5, eta = 0.8, fixed = "n1", difference = 0.2,
                  max = 300, batch_size = 100, seed = 1855, Bayes_pack = "bain")

set.seed(2916)
test_random <- round(runif(n = 5, min = 1, max = nrow(design_matrix_n2)))
results <- vector("list", 5)
index <- 1
for (i in test_random) {
    results[[index]] <- SSD_mult_CRT("homogeneity",
                                     effect_sizes = c(design_matrix_n2[i, 1], design_matrix_n2[i, 2]),
                                     n1 = design_matrix_n2[i, 10],
                                     n2 = 26,
                                     ndatasets = 100,
                                     out_specific_ICC = c(design_matrix_n2[i, 4], 0.1),
                                     intersubj_between_outICC = design_matrix_n2[i, 5],
                                     intrasubj_between_outICC = design_matrix_n2[i, 6],
                                     pmp_thresh = design_matrix_n2[i, 7],
                                     eta = design_matrix_n2[i, 8],
                                     fixed = as.character(design_matrix_n2[i, 9]),
                                     difference = design_matrix_n2[i, 3],
                                     max = 100,
                                     master.seed = 1629,
                                     Bayes_pack = "bain")
    index <- index + 1
}
# Data generation ---------------------------------------------------------------
source("multiv_data_generation.R")
source("helpers_functions.R")

options(error = recover)
a <- gen_multiv_data(30, 5, n2 = 16, c(0.3, 0.8), out_specific_ICC = 0.2,
                 intersubj_between_outICC = 0.022, intrasubj_between_outICC = 0.5,
                 n_outcomes = 2, seed = 40)

# Omnibus test------------------------------------------------------------------
set.seed(2916)
test_random <- round(runif(n = 5, min = 1, max = nrow(design_matrix_n2)))
results <- vector("list", 5)

## Serial
a <- SSD_mult_CRT("omnibus",
             effect_sizes = c(design_matrix_n2[i, 1], design_matrix_n2[i, 2]),
             n1 = design_matrix_n2[i, 10],
             n2 = 26,
             ndatasets = 100,
             out_specific_ICC = c(design_matrix_n2[i, 4], 0.1),
             intersubj_between_outICC = design_matrix_n2[i, 5],
             intrasubj_between_outICC = design_matrix_n2[i, 6],
             pmp_thresh = design_matrix_n2[i, 7],
             eta = design_matrix_n2[i, 8],
             fixed = as.character(design_matrix_n2[i, 9]),
             difference = design_matrix_n2[i, 3],
             max = 100,
             master.seed = 1629,
             Bayes_pack = "bain")

## Parallel
index <- 1
for (i in test_random) {
    results[[index]] <- SSD_mult_CRT("homogeneity",
                                     effect_sizes = c(design_matrix_n2[i, 1], design_matrix_n2[i, 2]),
                                     n1 = design_matrix_n2[i, 10],
                                     n2 = 26,
                                     ndatasets = 100,
                                     out_specific_ICC = c(design_matrix_n2[i, 4], 0.1),
                                     intersubj_between_outICC = design_matrix_n2[i, 5],
                                     intrasubj_between_outICC = design_matrix_n2[i, 6],
                                     pmp_thresh = design_matrix_n2[i, 7],
                                     eta = design_matrix_n2[i, 8],
                                     fixed = as.character(design_matrix_n2[i, 9]),
                                     difference = design_matrix_n2[i, 3],
                                     max = 100,
                                     master.seed = 1629,
                                     Bayes_pack = "bain")
    index <- index + 1
}

