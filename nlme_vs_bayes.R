################# Comparing nlme vs Bayesian estimation ####################

# packages
library(nlme)
library(brms)
library(ggmcmc)
library(MASS)
# Create dataset ------------------
n2 <- 40                                                   # number of clusters (should be even so that half of the clusters in control and other half in intervention)
n1 <- 5                                                    # common cluster size
cluster_id <- rep(1:n2, each = n1)                                  # ID for clusters
condit <- c(rep(0, each = (n1*n2)/2), rep(1,each = (n1*n2)/2))    # treatment condition

cov1 <- matrix(c(35, 15, 15, 40), 2, 2)                          # covariance matrix level 1
cov2 <- matrix(c(10, 4, 4, 12), 2, 2)

u <- mvrnorm(n2, mu = c(0,0), Sigma = cov2)
u1 <- u[,1]
u2 <- u[,2]

u1 <- rep(u1, each = n1)
u2 <- rep(u2, each = n1)

e <- mvrnorm(n1*n2, mu = c(0,0), Sigma = cov1)
e1 <- e[,1]
e2 <- e[,2]

y1 <- round(50 + 5*condit + u1 + e1)
y2 <- round(60 + 8*condit + u2 + e2)

dataset <- cbind(y1, y2, cluster_id, condit)
dataset <- as.data.frame(dataset)


## Long data format---------
dataset$obs_id <- 1:nrow(dataset)
long_data <- dataset %>%
    pivot_longer(
        cols = c(y1, y2),
        names_to = "outcome",
        values_to = "value"
    ) %>%
    mutate(
        outcome_fac = factor(outcome),
        outcome_num = as.numeric(outcome_fac)
    )
#------------------ Model with nlme ---------------------
model_nlme <- lme(
    value ~ 0 + outcome_fac + outcome_fac:condit,
    # Random intercept
    random = list(cluster_id = pdSymm(~ outcome_fac - 1)),
    # Correlation of observations within the clusters and within the same subject
    correlation = corSymm(form = ~ outcome_num | cluster_id/obs_id),
    # Variance of residuals per outcome
    weights = varIdent(form = ~ 1 | outcome_fac),
    data = long_data,
    control = lmeControl(opt = "optim", msMaxIter = 1000)
)
summary(model_nlme)

# variance-covariance of the fixed effects
var_cov_fixeff <- vcov(model_nlme)

# Variance-covariance matrix of random effects
var_random_eff <- getVarCov(model_nlme)

# Elements for ICCs
## Self-efficacy behaviour
var_u0_selfeff_b <- var_random_eff[1, 1]
var_e_selfeff_b <- unname(as.numeric(VarCorr(model_nlme)[3, 1]) * coef(model_nlme$modelStruct$varStruct, unconstrained = FALSE)**2)
total_var_selfeff_b <- var_u0_selfeff_b + var_e_selfeff_b
## Self-efficacy outcomes
var_u0_selfeff_o <- var_random_eff[2, 2]
var_e_selfeff_o <- as.numeric(VarCorr(model_nlme)[3, 1])
total_var_selfeff_o <- var_u0_selfeff_o + var_e_selfeff_o
# Both outcomes
cov_u0 <- var_random_eff[1, 2]
corr_e <- coef(model_nlme$modelStruct$corStruct, unconstrained = FALSE)
cov_e <- corr_e * sqrt(var_e_selfeff_b * var_e_selfeff_o)
### ICCs
# Self-efficacy behaviour
rho0_selfeff_o <- var_u0_selfeff_o / total_var_selfeff_o
# Self-efficacy behaviour
rho0_selfeff_b <- var_u0_selfeff_b / total_var_selfeff_b
rho1 <- cov_u0 / sqrt(total_var_selfeff_b * total_var_selfeff_o)
rho2 <- (cov_u0 + cov_e) / sqrt(total_var_selfeff_b * total_var_selfeff_o)
# Treatment effects
fix_eff <- fixed.effects(model_nlme)
treatm_eff_selfefficacy_b <- unname(fix_eff[3])/sqrt(total_var_selfeff_b)
treatm_eff_selfefficacy_o <- unname(fix_eff[4])/sqrt(total_var_selfeff_o)

#------------------- Model with brms ------------------------
# Bayes formula
bf_model <- bf(mvbind(y1, y2) ~ condit + (1|int|cluster_id)) + set_rescor(TRUE)

# Between the || is a label that identifies the random effect, if we include another
# random effect with the same label, then the package model them having a correlation
# between them.

fit_1 <- brm(bf_model, data = dataset, chains = 2, cores = 2)
summary(fit_1)

# Fixed effects
fixef(fit_1)
# Random effects sd-correlation matrix
summary(fit_1)$random
# Variance-covariance matrix of fixed effects
vcov(fit_1)


# 1. Variance-covariance of the fixed effects
# Returns the posterior covariance matrix of population-level effects
var_cov_fixeff_b <- vcov(fit_1)

# 2. Extract Variance-Covariance components
# VarCorr returns a list containing variances and correlations
vc <- VarCorr(fit_1)

## Random Effects (u0)
# Assuming 'id' is your grouping variable
var_u0_selfeff_b_bay <- vc$cluster_id$sd[1, "Estimate"]^2
var_u0_selfeff_o_bay <- vc$cluster_id$sd[2, "Estimate"]^2
cov_u0_bay           <- vc$cluster_id$cov[1, 1, 2]

## Residuals (e)
# In brms multivariate models, sigma is usually estimated per outcome
var_e_selfeff_b_bay  <- vc$residual$sd[1, "Estimate"]^2
var_e_selfeff_o_bay  <- vc$residual$sd[2, "Estimate"]^2
cov_e_bay            <- vc$residual$cov[1, 1, 2]

# 3. Calculate Totals
total_var_selfeff_b_bay <- var_u0_selfeff_b + var_e_selfeff_b
total_var_selfeff_o_bay <- var_u0_selfeff_o + var_e_selfeff_o

# 4. ICCs and Correlations
# ICC for Outcomes
rho0_selfeff_o_bay <- var_u0_selfeff_o / total_var_selfeff_o

# ICC for Behavior
rho0_selfeff_b_bay <- var_u0_selfeff_b / total_var_selfeff_b

# Between-group correlation
rho1_bay <- cov_u0 / sqrt(total_var_selfeff_b * total_var_selfeff_o)

# Total correlation (Combined)
rho2_bay <- (cov_u0 + cov_e) / sqrt(total_var_selfeff_b * total_var_selfeff_o)

# ========================== Conclusion =================================

#The two methods are almost identical. The differences between both are in  decimals.
# Specially, for the residual variance, where I found the largest difference between both methods
# of estimation.