################# EMPIRICAL APPLICATION WITH MOTIVATING EXAMPLE ################

# Libraries
library(dplyr)
library(nlme)
library(tidyr)
library(parallel)
library(nanoparquet)

# Required functions
source("find_n2_multiv.R")

# Data
self_affirmation <- read.csv(
  "~/Multivariate/exp_self-affirmation_data.csv",
    header = TRUE
)

#======================Data preparation=================================
head(self_affirmation)
str(self_affirmation)
summary(self_affirmation)
# From summary we can see that there are missing values marked as NA, but also as 99.
table(self_affirmation[ , "cw"]) # Number of clusters

## 99 to NA------
for (col in 1:ncol(self_affirmation)) {
    NA_rows <- which(self_affirmation[ , col] == 99)
    self_affirmation[ NA_rows, col] <- NA
}

## Create outcome variables----
### Self-worth positive index----
sw_items <- grep("^sw_", names(self_affirmation), value = TRUE)
sw_items <- sw_items[1:14]
# 1. Standard score
for (item in sw_items) {
    new_name <- paste0("z_", item)
    self_affirmation[[new_name]] <- scale(self_affirmation[[item]])[ , 1]
}

# 2. Average for each index
## positive
v_sum_pos <- vector()
for (i in 1:7) {
    v_sum_pos <- c(v_sum_pos, paste0("z_sw_0",i))
}
total_pos <- rowSums(self_affirmation[ , v_sum_pos], na.rm = TRUE)
self_affirmation$sw_positive <- total_pos/7
sd_sw_positive <- sd(self_affirmation$sw_positive, na.rm = TRUE)

## Negative
v_sum_neg <- vector()
for (i in 8:14) {
    if (i < 10) {
        v_sum_neg <- c(v_sum_neg, paste0("sw_0",i))
    } else {
        v_sum_neg <- c(v_sum_neg, paste0("sw_",i))
    }
}
total_neg <- rowSums(self_affirmation[ , v_sum_neg], na.rm = TRUE)
self_affirmation$sw_negative <- total_neg/7
sd_sw_negative <- sd(self_affirmation$sw_negative, na.rm = TRUE)

### Self-efficacy----
# Related to behaviours
v_sum <- vector()
for (i in 1:9) {
    v_sum <- c(v_sum, paste0("eff_beh_",i))
}
total_b <- rowSums(self_affirmation[ , v_sum], na.rm = TRUE)
self_affirmation$selfefficacy_b <- total_b/length(v_sum)
sd_self_efficacy_b <- sd(self_affirmation$selfefficacy_b, na.rm = TRUE)

# Related to outcomes
v_sum <- vector()
for (i in 1:7) {
    v_sum <- c(v_sum, paste0("eff_out_",i))
}
total_o <- rowSums(self_affirmation[ , v_sum], na.rm = TRUE)
self_affirmation$selfefficacy_o <- total_o/length(v_sum)
sd_self_efficacy_o <- sd(self_affirmation$selfefficacy_o, na.rm = TRUE)

## Long data format---------
self_affirmation$obs_id <- 1:nrow(self_affirmation)
long_data <- self_affirmation %>%
    pivot_longer(
        cols = c(selfefficacy_o, selfefficacy_b),
        names_to = "outcome",
        values_to = "value"
    ) %>%
    mutate(
        outcome_fac = factor(outcome),
        outcome_num = as.numeric(outcome_fac)
    )

#============ Fit bivariate model using nlme ===================

model_self_efficacy <- lme(
    value ~ 0 + outcome_fac + outcome_fac:treatment,
    # Random intercept
    random = list(cw = pdSymm(~ outcome_fac - 1)),
    # Correlation of observations within the clusters and within the same subject
    correlation = corSymm(form = ~ outcome_num | cw/obs_id),
    # Variance of residuals per outcome
    weights = varIdent(form = ~ 1 | outcome_fac),
    data = long_data,
    control = lmeControl(opt = "optim", msMaxIter = 1000)
)
summary(model_self_efficacy)

# variance-covariance of the fixed effects
var_cov_fixeff <- vcov(model_self_efficacy)

# Variance-covariance matrix of random effects
var_random_eff <- getVarCov(model_self_efficacy)

# Elements for ICCs
## Self-efficacy behaviour
var_u0_selfeff_b <- var_random_eff[1, 1]
var_e_selfeff_b <- unname(as.numeric(VarCorr(model_self_efficacy)[3, 1]) * coef(model_self_efficacy$modelStruct$varStruct, unconstrained = FALSE)**2)
total_var_selfeff_b <- var_u0_selfeff_b + var_e_selfeff_b
## Self-efficacy outcomes
var_u0_selfeff_o <- var_random_eff[2, 2]
var_e_selfeff_o <- as.numeric(VarCorr(model_self_efficacy)[3, 1])
total_var_selfeff_o <- var_u0_selfeff_o + var_e_selfeff_o
# Both outcomes
cov_u0 <- var_random_eff[1, 2]
corr_e <- coef(model_self_efficacy$modelStruct$corStruct, unconstrained = FALSE)
cov_e <- corr_e * sqrt(var_e_selfeff_b * var_e_selfeff_o)
### ICCs
# Self-efficacy behaviour
rho0_selfeff_o <- var_u0_selfeff_o / total_var_selfeff_o
# Self-efficacy behaviour
rho0_selfeff_b <- var_u0_selfeff_b / total_var_selfeff_b
rho1 <- cov_u0 / sqrt(total_var_selfeff_b * total_var_selfeff_o)
rho2 <- (cov_u0 + cov_e) / sqrt(total_var_selfeff_b * total_var_selfeff_o)
# Treatment effects
fix_eff <- fixed.effects(model_self_efficacy)
treatm_eff_selfefficacy_b <- unname(fix_eff[3])/sqrt(total_var_selfeff_b)
treatm_eff_selfefficacy_o <- unname(fix_eff[4])/sqrt(total_var_selfeff_o)

#======================== Sample Size Determination ============================
# Initial values----
pmp_thres <- c(0.9)
eta <- c(0.8, 0.9)
eff_sizes <- c(0.2, 0.3, 0.4, 0.5)
r <- 1
n1 <- 12 # Average cluster size
effects <- matrix(NA, ncol = 2, nrow = 16)
for (i in 1:length(eff_sizes)) {
    for (j in 1:length(eff_sizes)) {
        effects[r, ] <- c(eff_sizes[i], eff_sizes[j])
        r <- r + 1
    }
}
tests <- c("intersection-union", "omnibus", "homogeneity")
simulation_grid <- expand.grid(pmp_thres, eta, tests)
simulation_grid <- merge(as.data.frame(simulation_grid), as.data.frame(effects), by = NULL)
names(simulation_grid) <- c("pmp_thres", "eta", "test", "eff_1", "eff_2")
simulation_grid <- mutate(simulation_grid, seed = as.integer(sample(2^32/2, n())))
write_parquet(simulation_grid, "sim_application")
# 1. Intersection-union test and omnibus test
# H1 : β1SocBel > 0 & β1Self Ef f > 0
# H2 : β1SocBel > 0 & β1Self Ef f ≤ 0
# H3 : β1SocBel ≤ 0 & β1Self Ef f > 0
# H4 : β1SocBel ≤ 0 & β1Self Ef f ≤ 0

# 3. Homogeneity of treatment effects
# H1 : |β1SelfWorth − β1SelfEff | < Δ
# H2 : |β1SelfWorth − β1SelfEff | ≥ Δ

folder_results <- "application3"
if (!dir.exists(folder_results)) {
  dir.create(folder_results)
}

results <- mclapply(1:nrow(simulation_grid), function(Row) {
  
  
  # Start time
  start_time <- Sys.time()
  
  ssd_result <- SSD_mult_CRT(
    test = simulation_grid[Row, "test"],
    effect_sizes = c(simulation_grid[Row, "eff_1"],simulation_grid[Row, "eff_2"]),
    n1 = n1,
    n2 = 16,
    ndatasets = 500,
    out_specific_ICCs = c(rho0_selfeff_b, rho0_selfeff_o),
    intersubj_between_outICC = rho1,
    intrasubj_between_outICC = rho2,
    pmp_thresh = simulation_grid[Row, "pmp_thres"],
    eta = simulation_grid[Row, "eta"],
    fixed = "n1",
    max = 300,
    master.seed = simulation_grid[Row, "seed"],
    Bayes_pack = "bain", 
    difference = 0.3
  )
  
  # Save result
  end_time <- Sys.time()
  
  file_name <- file.path(folder_results, paste0("application", Row, ".RDS"))
  saveRDS(ssd_result, file = file_name)
  
  # Save running time
  running_time <- as.numeric(difftime(end_time,start_time, units = "mins"))
  time_name <- file.path(folder_results, paste0("time_", Row, ".RDS"))
  saveRDS(running_time, file = time_name)
  
  # Clean
  rm(ssd_results)
  gc()
  return(NULL)
}, mc.cores = 12)


#======================== Sample Size Determination with different ICCs ============================
# Initial values----
pmp_thres <- c(0.9, 0.95)
eta <- c(0.8, 0.9)
tests <- c("intersection-union", "omnibus", "homogeneity")
n1 <- 12 # Average cluster size
rho0 <- c(0.01, 0.05, 0.1)
pairs_rho0 <- matrix(NA, 9, 2)
r <- 1
for (i in 1:length(rho0)) {
  for (j in 1:length(rho0)) {
    pairs_rho0[r, ] <- c(rho0[i], rho0[j])
    r <- r + 1
  }
}
rho1 <- c(0.005, 0.025)
rho2 <- c(0.2, 0.5, 0.9)
simulation_grid <- expand.grid(pmp_thres, eta, tests, rho1, rho2)
simulation_grid <- merge(as.data.frame(simulation_grid), as.data.frame(pairs_rho0), by = NULL)
names(simulation_grid) <- c("pmp_thres", "eta", "test", "rho1", "rho2", "rho0.1", "rho0.2")
simulation_grid <- mutate(simulation_grid, seed = as.integer(sample(2^32/2, n())))


# 1. Intersection-union test and omnibus test
# H1 : β1SocBel > 0 & β1Self Ef f > 0
# H2 : β1SocBel > 0 & β1Self Ef f ≤ 0
# H3 : β1SocBel ≤ 0 & β1Self Ef f > 0
# H4 : β1SocBel ≤ 0 & β1Self Ef f ≤ 0

# 3. Homogeneity of treatment effects
# H1 : |β1SelfWorth − β1SelfEff | < Δ
# H2 : |β1SelfWorth − β1SelfEff | ≥ Δ

folder_results <- "application2"
if (!dir.exists(folder_results)) {
  dir.create(folder_results)
}

cores <- round(detectCores()*0.8)
results <- mclapply(1:nrow(simulation_grid), function(Row) {
  
  if (simulation_grid[Row, "test"] == "intersection-union") {
    name_ <- "IU_"
  } else if (simulation_grid[Row, "test"] == "omnibus") {
    name_ <- "omni_"
  } else if (simulation_grid[Row, "test"] == "homogeneity") {
    name_ <- "homog_"
  }
  
  
  # Start time
  start_time <- Sys.time()
  
  ssd_result <- SSD_mult_CRT(
    test = simulation_grid[Row, "test"],
    effect_sizes = c(treatm_eff_selfefficacy_b, (-1*treatm_eff_selfefficacy_o)),
    n1 = n1,
    n2 = 16,
    ndatasets = 500,
    out_specific_ICCs = c(simulation_grid[Row, "rho0.1"], simulation_grid[Row, "rho0.2"]),
    intersubj_between_outICC = simulation_grid[Row, "rho1"],
    intrasubj_between_outICC = simulation_grid[Row, "rho2"],
    pmp_thresh = simulation_grid[Row, "pmp_thres"],
    eta = simulation_grid[Row, "eta"],
    fixed = "n1",
    max = 300,
    master.seed = simulation_grid[Row, "seed"],
    Bayes_pack = "bain", 
    difference = 0.3
  )
  
  # Save result
  end_time <- Sys.time()
  
  file_name <- file.path(folder_results, paste0(name_, Row, ".RDS"))
  saveRDS(ssd_result, file = file_name)
  
  # Save running time
  running_time <- as.numeric(difftime(end_time,start_time, units = "mins"))
  time_name <- file.path(folder_results, paste0("time_", name_, Row, ".RDS"))
  saveRDS(running_time, file = time_name)
  
  # Clean
  rm(ssd_results)
  gc()
  return(NULL)
}, mc.cores = cores)

#---------------------------------------
Row <- 9
# Start time
start_time <- Sys.time()

ssd_result <- SSD_mult_CRT(
  test = simulation_grid[Row, "test"],
  effect_sizes = c(treatm_eff_selfefficacy_b, (-1*treatm_eff_selfefficacy_o)),
  n1 = n1,
  n2 = 16,
  ndatasets = 500,
  out_specific_ICCs = c( rho0_selfeff_b, rho0_selfeff_o),
  intersubj_between_outICC = rho1,
  intrasubj_between_outICC = rho2,
  pmp_thresh = simulation_grid[Row, "pmp_thres"],
  eta = simulation_grid[Row, "eta"],
  fixed = "n1",
  max = 300,
  master.seed = simulation_grid[Row, "seed"],
  Bayes_pack = "bain", 
  difference = 0.3
)

# Save result
end_time <- Sys.time()

file_name <- file.path(folder_results, paste0(name_, Row, ".RDS"))
saveRDS(ssd_result, file = file_name)

# Save running time
running_time <- as.numeric(difftime(end_time,start_time, units = "mins"))
time_name <- file.path(folder_results, paste0("time_", name_, Row, ".RDS"))
saveRDS(running_time, file = time_name)

# Clean
rm(ssd_results)
gc()

# Collect results ==============================================================
omnibus_index <- which(simulation_grid$test == "omnibus")
homogeneity_index <- which(simulation_grid$test == "homogeneity")
iu_index <- which(simulation_grid$test == "intersection-union")
application_omni <- collect_results(design_matrix = simulation_grid,
                               results_folder = folder_results, finding = "N2",
                               name_results = "application", test = "omnibus",
                               file_name = "application_omni", rows = omnibus_index,
                               save = FALSE)
application_homog <- collect_results(design_matrix = simulation_grid,
                                     results_folder = folder_results, finding = "N2",
                                     name_results = "application", test = "homogeneity",
                                     file_name = "application_homog", rows = homogeneity_index,
                                     save = FALSE)
application_iu <- collect_results(design_matrix = simulation_grid,
                                  results_folder = folder_results, finding = "N2",
                                  name_results = "application", test = "intersection-union",
                                  file_name = "application_iu", rows = iu_index, 
                                  save = FALSE)

run_again <- missing_rows("application3", name_pattern = "application", check_numbers = 1:96,
             underscore = FALSE)


results_ <- setdiff(homogeneity_index, run_again)
application_homog <- collect_results(design_matrix = simulation_grid,
                                     results_folder = folder_results, finding = "N2",
                                     name_results = "application", test = "homogeneity",
                                     file_name = "application_homog", rows = results_,
                                     save = FALSE)

application_omni <- application_omni[which(complete.cases(application_omni == TRUE)), ]
application_homog <- application_homog[which(complete.cases(application_homog == TRUE)), ]
application_iu <- application_iu[which(complete.cases(application_iu == TRUE)), ]

application_omni[which(application_omni$eff_1 == 0.2 & (application_omni$eff_2 == 0.2 | application_omni$eff_2 == 0.3)), ]
application_homog[which(application_homog$eff_1 == 0.2 & (application_homog$eff_2 == 0.2 | application_homog$eff_2 == 0.3)), ]
application_iu[which(application_iu$eff_1 == 0.2 & (application_iu$eff_2 == 0.2 | application_iu$eff_2 == 0.3)), ]

