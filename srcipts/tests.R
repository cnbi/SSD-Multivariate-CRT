############# TESTS###########################################

source("find_n2_multiv.R")

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

# Hypothesis true
start_time <- Sys.time()
SSD_mult_CRT("intersection-union",
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
             Bayes_pack = "BFpack")
end_time <- Sys.time()
as.numeric(difftime(end_time, start_time, units = "mins"))

start_time <- Sys.time()
SSD_mult_CRT("intersection-union",
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
             Bayes_pack = "BFpack")

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