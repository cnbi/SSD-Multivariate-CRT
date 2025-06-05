# Data generation

# ndatasets: Number of datasets to generate
# n1: Cluster size. We assume equal size for all clusters
# n2: Number of clusters. It is the total number of clusters, considering that we have two treatment conditions, the number must be even.
# effect_sizes: Vector with the effect sizes.
# intersubject_iccs: A matrix with the intersubject intraclass correlations. The diagonal corresponds to the endpoint-specific icc and outside the diagonal corresponds to between-endpoint icc.
# intrasubject_icc: A matrix with the intrasubject between-endpoint intraclass correlation outside the diagonal.The diagonal contains 1.
# var_y: Vector with the variances of each outcome variable.

# Libraries --------------------------------------------------------
library(nlme) #MLMM
library(lme4)
library(MASS) #Multivariate data generation
library(reshape) #long format

source("EM_algorithm.R")

# Data generation --------------------------------------------------
#Starting values
ndatasets <- 10
n1 <- 30
n2 <- 100
effect_sizes <- c(0.5, 0.3)
n_outcomes <- length(effect_sizes)
icc <- 0.05

# Matrix with intersubject intraclass correlation coefficients
intersubj_iccs <- matrix(NA, n_outcomes, n_outcomes)
diag_endpoint_icc <- c(icc, icc) #rho_0 of the two outcomes
intersubj_iccs <- diag(diag_endpoint_icc)
values_intersubj_mat <- c(0.005) #rho_1: order y2y1, y3y1, y3y2
i <- 0
for (column in 1:(n_outcomes - 1)) {
    for (row in 2:n_outcomes) {
        if (column != row) {
            i <- i + 1
            intersubj_iccs[row, column] <- values_intersubj_mat[i] #rho_1
            intersubj_iccs[column, row] <- values_intersubj_mat[i] #rho_1
        }
    }
}
y_names <- c("y1", "y2")
colnames(intersubj_iccs) <- y_names
rownames(intersubj_iccs) <- y_names

# Matrix with intrasubject between-endpoints intraclass correlation coefficients
intrasubj_icc <- matrix(NA, n_outcomes, n_outcomes)
intrasubj_icc <- diag(1, n_outcomes)
values_intrasubj_mat <- rep(0.2, (n_outcomes*(n_outcomes - 1)/2)) #rho_1: order y2y1, y3y1, y3y2
i <- 0
for (column in 1:(n_outcomes - 1)) {
    for (row in 2:n_outcomes) {
        if (column != row) {
            i <- i + 1
            intrasubj_icc[row, column] <- values_intrasubj_mat[i] #rho_1
            intrasubj_icc[column, row] <- values_intrasubj_mat[i] #rho_1
        }
    }
}
colnames(intrasubj_icc) <- y_names
rownames(intrasubj_icc) <- y_names
var_y <- rep(1, n_outcomes) # Marginal variance=total variance

# Objects to save results
output_multilevel <- vector(mode = "list", length = ndatasets)
data_list <- vector(mode = "list", length = ndatasets)

# Data generation
## Treatment condition
id_cluster <- rep(1:n2, each = n1)
id_subj <- seq(1:(n1 * n2))
condition <- rep(c(0, 1), each = n1 * n2 / 2)

## Covariance matrices specification 
#Sigma e
sigma_e <- diag(1 - diag(intersubj_iccs) * var_y)
i <- 0
for (column in 1:(n_outcomes - 1)) {
    for (row in 2:n_outcomes) {
        if (column != row) {
            i <- i + 1
            sigma_e[row, column] <- sqrt(var_y[row]) * sqrt(var_y[column]) * (intrasubj_icc[row, column] - intersubj_iccs[row, column])
            sigma_e[column, row] <- sqrt(var_y[row]) * sqrt(var_y[column]) * (intrasubj_icc[column, row] - intersubj_iccs[column, row])
        }
    }
}

## Sigma u0
sigma_u0 <- diag(diag(intersubj_iccs) * var_y)
i <- 0
for (column in 1:(n_outcomes - 1)) {
    for (row in 2:n_outcomes) {
        if (column != row) {
            i <- i + 1
            sigma_u0[row, column] <- sqrt(var_y[row]) * sqrt(var_y[column]) * intersubj_iccs[row, column]
            sigma_u0[column, row] <- sqrt(var_y[row]) * sqrt(var_y[column]) * intersubj_iccs[column, row]
        }
    }
}

# Random effects
set.seed(425)
e <- MASS::mvrnorm(n1*n2, rep(0, n_outcomes), sigma_e)
u0 <- MASS::mvrnorm(n2, rep(0, n_outcomes), sigma_u0)
u0 <- do.call(rbind, replicate(n1, u0, simplify = FALSE))


# betay1 <- effect_sizes[1] * sqrt(sigma_u0[1, 1])
# betay2 <- effect_sizes[2] * sqrt(sigma_u0[2, 2])
# y1.b <- betay1 * condition + u0[, 1] + e[, 1]
# y2.b <- betay2 * condition + u0[, 2] + e[, 2]
# df <- cbind(id_subj, id_cluster, condition, y1.b, y2.b)

# Means
# Dummy variables for no intercept model
intervention <- condition
control <- 1 - intervention
mean_control <- 0
mean_interv <- effect_sizes
y1. <- mean_control * control + mean_interv[1] * intervention + u0[, 1] + e[, 1]
y2. <- mean_control * control + mean_interv[2] * intervention + u0[, 2] + e[, 2]

# my_data <- cbind(id_subj, id_cluster, intervention, control, y1,y1.b, y2, y2.b, condition)
my_data. <- cbind(id_subj, id_cluster, intervention, control, y1., y2., condition)

# Multivariate multilevel model-------------------------------------------------
library(nlme)
library(mvtnorm)
library(numDeriv)

source("EM_algorithm.R")
source("loglikelihood.R")

colnames(my_data.)[2] <- "cluster"
estimations. <- EM.estim2(as.data.frame(my_data.), as.formula('y1. ~ condition'), as.formula('y2. ~ condition'), 
          maxiter = 500, verbose = TRUE)


# Extract the elements for outcome 1
estimates_y1 <- estimations$theta$zeta[1:2]
names(estimates_y1) <- c("intercept", "slope")

# Extract the elements for outcome 2
estimates_y2 <- estimations$theta$zeta[3:4]
names(estimates_y2) <- c("Intercept", "slope")

# Bayes factor -------------------------------------------------

## Y1
# Complexities
comp1 <- .5
comp2 <- .5

# Fit
fit2 <- pnorm(0, mean = estimates_y1[2], sd = sqrt(estimations$var_cov[2,2]))
fit1 <- 1 - fit2

# Calculation BFs
# H1: control<treatment (slope>0)
# H2: control>treatment (slope<0)
AAFBF1u <- (fit1/comp1) 
AAFBF2u <- (fit2/comp2) 
AAFBF12 <- AAFBF1u/AAFBF2u 
AAFBF21 <- 1/AAFBF12

#Calculation of PMPs
pmp1 <- AAFBF1u/(AAFBF1u + AAFBF2u)
pmp2 <- 1 - pmp1
output <- list(bf.12 = AAFBF12, bf.21 = AAFBF21, pmp1 = pmp1, pmp2 = pmp2)

# Using bain
library(bain)
effective_n <- (n1*n2) / (1 + (n1 - 1) * icc) / 2
bain(estimates_y1, "slope<0; slope>0", n = effective_n,
     Sigma = estimations$hessian_method[1:2, 1:2]) # The result is similar to my code

# Using BFpack
library("BFpack")

bf_y1 <- BF(estimates_y1, Sigma = estimations$var_cov[1:2, 1:2], n = effective_n, 
            hypothesis = "slope<0; slope>0")
# Because I am working with hypotheses with only inequality constraints, I have compare
# my results with  the 6th column (BF>). This column show the BF of the inequality 
# constained hypothesis against the unconstrained hypothesis. 


## Y2
# Fit
fit2 <- pnorm(0, mean = estimates_y2[2], sd = sqrt(estimations$var_cov[4,4]))
fit1 <- 1 - fit2

# Calculation BFs
# H1: control<treatment
# H2: control>treatment
AAFBF1u <- (fit1/comp1)
AAFBF2u <- (fit2/comp2)
AAFBF12 <- AAFBF1u/AAFBF2u
AAFBF21 <- 1/AAFBF12

#Calculation of PMPs
pmp1 <- AAFBF1u/(AAFBF1u + AAFBF2u)
pmp2 <- 1 - pmp1
output <- list(bf.12 = AAFBF12, bf.21 = AAFBF21, pmp1 = pmp1, pmp2 = pmp2)

# Using bain
library(bain)

bain(estimates_y2, "slope<0; slope>0", n = effective_n,
     Sigma = estimations$hessian_method[3:4, 3:4]) # The result is similar to my code

# Using BFpack
library("BFpack")

bf_y2 <- BF(estimates_y2, Sigma = estimations$var_cov[3:4, 3:4], n = effective_n, 
            hypothesis = "slope<0; slope>0")
# Because I am working with hypotheses with only inequality constraints, I have compare
# my results with  the 6th column (BF>). This column show the BF of the inequality 
# constrained hypothesis against the unconstrained hypothesis. 


bf_y2hess <- BF(estimates_y2, Sigma = estimations$hessian_method[3:4, 3:4], n = effective_n, 
            hypothesis = "slope<0; slope>0")

# Homogeneity Test -------------------------------------------
# $H_0$: Effect sizes (slopes) or treatment effects are equal (across endpoints). 
# $H_1$: At least one effect size is different compared with the other effect sizes.
## Using bain
estimations_y1y2 <- c(estimates_y1, estimates_y2)
names(estimations_y1y2) <- c("Intercept1", "Slope1", "Intercept2", "Slope2") 
bain(estimations_y1y2, "Slope1=Slope2", n = effective_n,
     Sigma = estimations$hessian_method)

## Uisng BFpack
bf_y2hess <- BF(estimations_y1y2, Sigma = estimations$hessian_method, n = effective_n, 
                hypothesis = "Slope1=Slope2")
summary(bf_y2hess)


# Intersection--union test ------------------------------------
# $H_0$: At least one or both effect sizes are zero.
# $H_1$: Effect sizes are larger than zero.
# This means that now we have to test that the effect sizes are larger than zero and the Bayes factor indicating this should be at least the Bayes factor threshold.

intersection_union <- function(estimations_y1y2, covariances_matrix, effective_n, Bayes_pack) {
    result <- FALSE
    if(Bayes_pack == "bain"){
        Bf <- bain(estimations_y1y2, "Slope1>0; Slope2>0", n = effective_n,
                   Sigma = covariances_matrix)
    } else if (Bayes_pack == "BFpack"){
        Bf <- BF(estimations_y1y2, Sigma = estimations$hessian_method, n = effective_n, 
                 hypothesis = "Slope1>0; Slope2>0")
    }
    browser()
    #ifelse(Bf$)
}