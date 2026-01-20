################# EMPIRICAL APPLICATION WITH MOTIVATING EXAMPLE ################

# Data
self_affirmation <- read.csv(
    "C:\\Users\\barra006\\Nextcloud\\Second_project\\self-affirmation\\exp_self-affirmation_data.csv",
    header = TRUE
)
head(self_affirmation)
str(self_affirmation)
summary(self_affirmation)
# From summary we can see that there are missing values marked as NA, but also as 99.
table(self_affirmation[ , "cw"]) # Number of clusters
# Required functions
source("find_n2_multiv.R")

# Initial values
intersubj_icc <- c(0.005, 0.025)
intrasubj_icc <- c(0.2, 0.5)
pmp_thres <- c(0.9, 0.95)
eta <- c(0.8, 0.9)
simulation_grid <- expand.grid(intersubj_icc, intrasubj_icc, pmp_thres, eta)
names(simulation_grid) <- c("rho_1", "rho_2", "pmp_thres", "eta")

# 1. Intersection-union test ####
# According to the section 5. Test for Multiple Outcomes in the test, the hypotheses
# to test are:
# H1 : β1SocBel > 0 & β1Self Ef f > 0
# H2 : β1SocBel > 0 & β1Self Ef f ≤ 0
# H3 : β1SocBel ≤ 0 & β1Self Ef f > 0
# H4 : β1SocBel ≤ 0 & β1Self Ef f ≤ 0

for (Row in 1:nrow(simulation_grid)) {
    number_of_clusters[[Row]] <- SSD_mult_CRT(
        test = "intersection-union",
        effect_sizes = c(0.49, 0.49),
        n1 = 20,
        n2 = 15,
        ndatasets = 500,
        out_specific_ICCs = 0.05,
        intersubj_between_outICC = simulation_grid[Row, "rho_1"],
        intrasubj_between_outICC = simulation_grid[Row, "rho_2"],
        pmp_thresh = simulation_grid[Row, "pmp_thres"],
        eta = simulation_grid[Row, "eta"],
        fixed = "n1",
        max = 300,
        master.seed = 14012026,
        Bayes_pack = "bain"
    )
}


# 2. Omnibus test ####
# H1 : β1SelfWorth > 0&β1SelfEff ≤ 0
# H2 : β1SelfWorth ≤ 0&β1SelfEff > 0
# H3 : β1SelfWorth > 0&β1SelfEff > 0
# H4 : β1SelfWorth ≤ 0&β1SelfEff ≤ 0
for (Row in 1:nrow(simulation_grid)) {
    number_of_clusters[[Row]] <- SSD_mult_CRT(
        test = "omnibus",
        effect_sizes = c(0.49, 0.49),
        n1 = 20,
        n2 = 15,
        ndatasets = 500,
        out_specific_ICCs = 0.05,
        intersubj_between_outICC = simulation_grid[Row, "rho_1"],
        intrasubj_between_outICC = simulation_grid[Row, "rho_2"],
        pmp_thresh = simulation_grid[Row, "pmp_thres"],
        eta = simulation_grid[Row, "eta"],
        fixed = "n1",
        max = 300,
        master.seed = 14012026,
        Bayes_pack = "bain"
    )
}


# 3. Homogeneity of treatment effects #####
# H1 : |β1SelfWorth − β1SelfEff | < Δ
# H2 : |β1SelfWorth − β1SelfEff | ≥ Δ
for (Row in 1:nrow(simulation_grid)) {
    number_of_clusters[[Row]] <- SSD_mult_CRT(
        test = "omnibus",
        effect_sizes = c(0.49, 0.49),
        n1 = 20,
        n2 = 15,
        ndatasets = 500,
        out_specific_ICCs = 0.05,
        intersubj_between_outICC = simulation_grid[Row, "rho_1"],
        intrasubj_between_outICC = simulation_grid[Row, "rho_2"],
        pmp_thresh = simulation_grid[Row, "pmp_thres"],
        eta = simulation_grid[Row, "eta"],
        fixed = "n1",
        max = 300,
        master.seed = 14012026,
        Bayes_pack = "bain",
        difference = 0.3
    )
}

###############################################
# Sample Size Determination for a Replication #
###############################################

# Create outcome variables
## Self-worth positive index
v_sum <- vector()
for (i in 1:7) {
    v_sum <- c(v_sum, paste0("sw_0",i))
}
self_affirmation$sw_positive <- rowSums(self_affirmation[ , v_sum], na.rm = TRUE)
v_sum <- vector()
for (i in 8:14) {
    if (i < 10) {
        v_sum <- c(v_sum, paste0("sw_0",i))
    } else {
        v_sum <- c(v_sum, paste0("sw_",i))
    }
}
sd_sw_positive <- sd(self_affirmation$sw_positive, na.rm = TRUE)

## Self-worth negative index
self_affirmation$sw_negative <- rowSums(self_affirmation[ , v_sum], na.rm = TRUE)

## Stress
self_affirmation$stress <- rowSums(self_affirmation[ , c("exp_2", "exp_4", "exp_6")], na.rm = TRUE)

## Job search behavioural self-efficacy
v_sum <- vector()
for (i in 1:9) {
    v_sum <- c(v_sum, paste0("eff_beh_",i))
}
self_affirmation$job_behav_selfefficacy <- rowSums(self_affirmation[ , v_sum], na.rm = TRUE)
sd_self_efficacy <- sd(self_affirmation$job_behav_selfefficacy, na.rm = TRUE)

# Societal Belonging
sd_soc_bel <- sd(self_affirmation$soc_belonging, na.rm = TRUE) 

# Initial values
intersubj_icc <- c(0.005, 0.025)
intrasubj_icc <- c(0.2, 0.5)
pmp_thres <- c(0.9, 0.95)
eta <- c(0.8, 0.9)
simulation_grid <- expand.grid(intersubj_icc, intrasubj_icc, pmp_thres, eta)
names(simulation_grid) <- c("rho_1", "rho_2", "pmp_thres", "eta")

# Transform from OLS estimation to Cohen's d
d_soc_bel <- -0.448/sd_soc_bel
d_self_eff <- 0.015/sd_self_efficacy
d_sw_positive <- -0.412/sd_sw_positive

# 1. Intersection-union test
# H1 : β1SocBel > 0 & β1Self Ef f > 0
# H2 : β1SocBel > 0 & β1Self Ef f ≤ 0
# H3 : β1SocBel ≤ 0 & β1Self Ef f > 0
# H4 : β1SocBel ≤ 0 & β1Self Ef f ≤ 0
for (Row in 1:nrow(simulation_grid)) {
    number_of_clusters[[Row]] <- SSD_mult_CRT(
        test = "intersection-union",
        effect_sizes = c(d_soc_bel, d_self_eff),
        n1 = 20,
        n2 = 15,
        ndatasets = 500,
        out_specific_ICCs = 0.05,
        intersubj_between_outICC = simulation_grid[Row, "rho_1"],
        intrasubj_between_outICC = simulation_grid[Row, "rho_2"],
        pmp_thresh = simulation_grid[Row, "pmp_thres"],
        eta = simulation_grid[Row, "eta"],
        fixed = "n1",
        max = 300,
        master.seed = 14012026,
        Bayes_pack = "bain"
    )
}

# 2. Omnibus test
# H1 : β1SelfWorth > 0&β1SelfEff ≤ 0
# H2 : β1SelfWorth ≤ 0&β1SelfEff > 0
# H3 : β1SelfWorth > 0&β1SelfEff > 0
# H4 : β1SelfWorth ≤ 0&β1SelfEff ≤ 0
for (Row in 1:nrow(simulation_grid)) {
    number_of_clusters[[Row]] <- SSD_mult_CRT(
        test = "omnibus",
        effect_sizes = c(d_sw_positive, d_self_eff),
        n1 = 20,
        n2 = 15,
        ndatasets = 500,
        out_specific_ICCs = 0.05,
        intersubj_between_outICC = simulation_grid[Row, "rho_1"],
        intrasubj_between_outICC = simulation_grid[Row, "rho_2"],
        pmp_thresh = simulation_grid[Row, "pmp_thres"],
        eta = simulation_grid[Row, "eta"],
        fixed = "n1",
        max = 300,
        master.seed = 14012026,
        Bayes_pack = "bain"
    )
}

# 3. Homogeneity of treatment effects
# H1 : |β1SelfWorth − β1SelfEff | < Δ
# H2 : |β1SelfWorth − β1SelfEff | ≥ Δ
for (Row in 1:nrow(simulation_grid)) {
    number_of_clusters[[Row]] <- SSD_mult_CRT(
        test = "omnibus",
        effect_sizes = c(d_sw_positive, d_self_eff),
        n1 = 20,
        n2 = 15,
        ndatasets = 500,
        out_specific_ICCs = 0.05,
        intersubj_between_outICC = simulation_grid[Row, "rho_1"],
        intrasubj_between_outICC = simulation_grid[Row, "rho_2"],
        pmp_thresh = simulation_grid[Row, "pmp_thres"],
        eta = simulation_grid[Row, "eta"],
        fixed = "n1",
        max = 300,
        master.seed = 14012026,
        Bayes_pack = "bain",
        difference = 0.3
    )
}
