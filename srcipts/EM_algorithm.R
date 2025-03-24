
###################### EM algorithm #########3
# data
# formula: dep ~ ind
# maxinter
# epsilon


# function to perform EM estimation with K=2 outcomes
EM.estim2 <- function(data, formula1, formula2, maxiter = 500, epsilon = 1e-4, 
                      verbose = FALSE){
    # data: source data set
    # formula1: fit mixed model with outcome 1 and the treatment arm. e.g formula1 <-as.formula(  "out1 ~ arm")
    # formula2: fit mixed model with outcome 2 and the treatment arm. e.g formula1 <-as.formula(  "out2 ~ arm")
    
    # fit mixed model to initialize parameters
    fm1 <- lme(formula1, random = ~ 1|cluster, data = data) #y1~condition
    fm2 <- lme(formula2, random = ~ 1|cluster, data = data) #y2~condition
    zeta <- as.numeric(c(fm1$coefficients$fixed, fm2$coefficients$fixed)) #fixed effects coefficients
    beta1 <- zeta[1:2] 
    beta2 <- zeta[3:4]
    
    n1s <- as.numeric(table(data$cluster)) #Vector with n1 for every cluster
    s2phi1 <- VarCorr(fm1)[1, 1]  # Random intercepts model 1
    s2phi2 <- VarCorr(fm2)[1, 1]  # Random intercepts model 2
    SigmaPhi <- diag(c(s2phi1, s2phi2))  #Variance-covariance matrix for random intercepts
    InvS2Phi <- solve(SigmaPhi)
    
    s2e1 <- VarCorr(fm1)[2, 1] # Residual variance model 1
    s2e2 <- VarCorr(fm2)[2, 1] # Residual variance model 2
    SigmaE <- diag(c(s2e1, s2e2)) # Variance-covariance matrix for residuals
    InvS2E <- solve(SigmaE)
    
    # Outcomes and treatment arm names
    out1 <- all.vars(formula1)[1] 
    out2 <- all.vars(formula2)[1]
    arm <- all.vars(formula1)[2]

    Y <- as.matrix(data[,c(out1, out2)])            # Matrix with only outcomes
    ID <- as.numeric(data$cluster)                  # ID for each cluster
    n2 <- length(unique(ID))                        # n2
    facz <- as.factor(data[,arm])               # Treatment arm variable as factor
    z <- as.numeric(facz) - 1                   # Control 0 and treatment 1
    X <- as.matrix(cbind(1, z))                 # design matrix for fixed effects
    n_outcomes <- 2
    ESSphi1 <- matrix(0, n2, n_outcomes)        # Expected random effects
    ESSphi2 <- array(0, c(n_outcomes, n_outcomes, n2))  # Covariance of random effects
    
    delta <- 2*epsilon                           # Change in log-likelihood
    converge <- 0
    
    # log likelihood
    theta <- c(zeta, c(SigmaPhi[!lower.tri(SigmaPhi)]),c(SigmaE[!lower.tri(SigmaE)]))  # Parameters
    LLold <- loglik(theta, n2, Y, X, ID, n1s, outcomes = 2)  # Previous log-likelihood
    
    
    niter <- 1
    while ((niter <= maxiter) & (abs(delta) > epsilon)) {
        
        # Expectation step
        for (j in 1:n2) {
            Yj <- Y[ID == j,,drop = FALSE]
            Xj <- X[ID == j,,drop = FALSE]
            residj <- Yj - cbind(Xj %*% beta1, Xj %*% beta2)
            Vj <- solve(InvS2Phi + n1s[j]*InvS2E)
            Muj <- as.numeric(Vj %*% InvS2E %*% colSums(residj))
            Nujj <- Vj + tcrossprod(Muj)
            ESSphi1[j, ] <- Muj
            ESSphi2[, , j] <- Nujj
        }
        
        # Maximization step - phi
        SigmaPhi <- apply(ESSphi2, 1:2, sum)/n2
        InvS2Phi <- solve(SigmaPhi)
        
        # Maximization step - zeta
        # Simplify the expression analytically, and obtain simple expression!
        XXt <- crossprod(X)
        Vzeta <- solve(kronecker(InvS2E, XXt))
        rzeta1 <- t(X) %*% (Y[, 1] - ESSphi1[ID, 1])
        rzeta2 <- t(X) %*% (Y[, 2] - ESSphi1[ID, 2])
        zeta <- Vzeta %*% rbind(InvS2E[1, 1] * rzeta1 + InvS2E[1, 2] * rzeta2,
                                InvS2E[2, 1] * rzeta1 + InvS2E[2, 2] * rzeta2)
        zeta <- c(zeta)
        beta1 <- zeta[1:2]
        beta2 <- zeta[3:4]
        
        # Maximization step - epsilon
        re <- Y - cbind(X %*% beta1, X %*% beta2)
        rss <- crossprod(re) + rowSums(sweep(ESSphi2, 3, n1s, FUN = "*"),dims = 2) -
            crossprod(ESSphi1, rowsum(re, ID)) - crossprod(rowsum(re, ID), ESSphi1)
        SigmaE <- rss/sum(n1s)
        # SigmaE <- diag(diag(SigmaE))
        InvS2E <- solve(SigmaE)
        
        # whether the algorithm converges
        # theta = c(zeta,c(SigmaPhi[!lower.tri(SigmaPhi)]),diag(SigmaE))
        theta  <- c(zeta, c(SigmaPhi[!lower.tri(SigmaPhi)]), c(SigmaE[!lower.tri(SigmaE)]))
        LLnew <- loglik(theta, n2, Y, X, ID, n1s, 2)
        delta <- abs(LLnew - LLold)
        LLold <- LLnew
        converge <- (abs(delta) <= epsilon)
        niter <- niter + 1
        if (verbose) cat(paste('iter=',niter), '\t', 
                         paste('param.error=',epsilon), '\t',
                         paste('loglik=',LLnew),'\n');
        
    }
    #Variance-covariance matrix
    hessian_matrix <- hessian(neg_loglik, theta, n2 = n2, Y = Y, X = X, ID = ID, n1s = n1s, outcomes = 2)
    var_cov <- solve(hessian_matrix)
    param <- list(theta = list(zeta = zeta, SigmaE = SigmaE, SigmaPhi = SigmaPhi),
                  loglik = LLnew, eps = epsilon, iter = niter,
                  var_cov = var_cov)
    return(param)
}

# # function to perform EM estimation with K=3 outcomes
# EM.estim3 <- function(data, formula1, formula2, formula3, maxiter = 500, epsilon = 1e-4
#                       , verbose = FALSE){
#     # data: source data set
#     # formula1: fit mixed model with outcome 1 and the treatment arm. e.g formula1 <-as.formula(  "out1 ~ arm")
#     # formula2: fit mixed model with outcome 2 and the treatment arm. e.g formula1 <-as.formula(  "out2 ~ arm")
#     # formula2: fit mixed model with outcome 3 and the treatment arm. e.g formula1 <-as.formula(  "out3 ~ arm")
#     
#     #Libraries
#     library(nlme)
#     library(mvtnorm)
#     library(numDeriv)
#     
#     #Functions
#     source("loglikelihood.R")
#     
#     # fit mixed model to initialize parameters
#     fm1 <- lme(formula1, random = ~ 1|cluster, data = data)
#     fm2 <- lme(formula2, random = ~ 1|cluster, data = data)
#     fm3 <- lme(formula3, random = ~ 1|cluster, data = data)
#     
#     n_outcomes <- 3
#     zeta <- as.numeric(c(fm1$coefficients$fixed, fm2$coefficients$fixed, fm3$coefficients$fixed))
#     beta1 <- zeta[1:2]
#     beta2 <- zeta[3:4]
#     beta3 <- zeta[5:6]
#     
#     n1s <- as.numeric(table(data$cluster))
#     
#     s2phi1 <- VarCorr(fm1)[1, 1]
#     s2phi2 <- VarCorr(fm2)[1, 1]
#     s2phi3 <- VarCorr(fm3)[1, 1]
#     
#     SigmaPhi <- diag(c(s2phi1, s2phi2, s2phi3))
#     InvS2Phi <- solve(SigmaPhi)
#     
#     s2e1 <- VarCorr(fm1)[2, 1]
#     s2e2 <- VarCorr(fm2)[2, 1]
#     s2e3 <- VarCorr(fm3)[2, 1]
#     
#     SigmaE <- diag(c(s2e1, s2e2, s2e3))
#     InvS2E <- solve(SigmaE)
#     
#     out1 <- all.vars(formula1)[1]
#     out2 <- all.vars(formula2)[1]
#     out3 <- all.vars(formula3)[1]
#     
#     arm <- all.vars(formula1)[2]
#     
#     Y <- as.matrix(data[,c(out1, out2, out3)])
#     ID <- as.numeric(data$cluster)
#     n2 <- length(unique(ID))
#     #X <- as.matrix(cbind(1, data[,"arm"])) # design matrix
#     facz <- as.factor(data[,arm])
#     z <- as.numeric(facz) - 1
#     X <- as.matrix(cbind(1, z)) # design matrix
#     
#     ESSphi1 <- matrix(0, n2, n_outcomes)
#     ESSphi2 <- array(0, c(n_outcomes, n_outcomes, n2))
#     
#     
#     #maxiter=500
#     #epsilon=1e-4
#     delta <- 2 * epsilon
#     max_modi <- 20
#     
#     converge <- 0
#     
#     # log likelihood
#     theta <- c(zeta, c(SigmaPhi[!lower.tri(SigmaPhi)]),c(SigmaE[!lower.tri(SigmaE)]))
#     LLold <- loglik(theta, n2, Y, X, ID, n1s, n_outcomes = 3)
#     
#     
#     niter <- 1
#     while ((niter <= maxiter) & (abs(delta) > epsilon)) {
#         
#         # Expectation step
#         for (j in 1:n2) {
#             Yj <- Y[ID == j, , drop = FALSE]
#             Xj <- X[ID == j, , drop = FALSE]
#             residj <- Yj - cbind(Xj %*% beta1, Xj %*% beta2, Xj %*% beta3)
#             Vj <- solve(InvS2Phi + n1s[j]*InvS2E)
#             Muj <- as.numeric(Vj %*% InvS2E %*% colSums(residj))
#             Nujj <- Vj + tcrossprod(Muj)
#             ESSphi1[j, ] <- Muj
#             ESSphi2[, , j] <- Nujj
#         }
#         
#         # Maximization step - phi
#         SigmaPhi <- apply(ESSphi2, 1:2, sum)/n2
#         InvS2Phi <- solve(SigmaPhi)
#         
#         # Maximization step - zeta
#         # Simplify the expression analytically, and obtain simple expression!
#         XXt <- crossprod(X)
#         Vzeta <- solve(kronecker(InvS2E, XXt))
#         rzeta1 <- t(X) %*% (Y[, 1] - ESSphi1[ID, 1])
#         rzeta2 <- t(X) %*% (Y[, 2] - ESSphi1[ID, 2])
#         rzeta3 <- t(X) %*% (Y[, 3] - ESSphi1[ID, 3])
#         
#         zeta <- Vzeta %*% rbind(InvS2E[1, 1]*rzeta1 + InvS2E[1, 2]*rzeta2 + InvS2E[1, 3]*rzeta3,
#                                 InvS2E[2, 1]*rzeta1 + InvS2E[2, 2]*rzeta2 + InvS2E[2, 3]*rzeta3,
#                                 InvS2E[3, 1]*rzeta1 + InvS2E[3, 2]*rzeta2 + InvS2E[3, 3]*rzeta3
#         )
#         zeta <- c(zeta)
#         beta1 <- zeta[1:2]
#         beta2 <- zeta[3:4]
#         beta3 <- zeta[5:6]
#         
#         # Maximization step - epsilon
#         re <- Y - cbind(X %*% beta1, X %*% beta2, X %*% beta3)
#         rss <- crossprod(re) + rowSums(sweep(ESSphi2, 3, n1s, FUN = "*"), dims = 2) -
#             crossprod(ESSphi1, rowsum(re, ID)) - crossprod(rowsum(re, ID),ESSphi1)
#         SigmaE <- rss/sum(n1s)
#         # SigmaE <- diag(diag(SigmaE))
#         InvS2E <- solve(SigmaE)
#         
#         # whether the algorithm converges
#         # theta = c(zeta, c(SigmaPhi[!lower.tri(SigmaPhi)]),diag(SigmaE))
#         theta = c(zeta, c(SigmaPhi[!lower.tri(SigmaPhi)]),c(SigmaE[!lower.tri(SigmaE)]))
#         #theta = c(zeta, diag(SigmaPhi),SigmaPhi[1, 2], c(SigmaE[!lower.tri(SigmaE)]))
#         
#         LLnew <- loglik(theta, n2, Y, X, ID, n1s, outcomes = 3)
#         delta <- abs(LLnew - LLold)
#         LLold <- LLnew
#         converge <- (abs(delta) <= epsilon)
#         niter <- niter + 1
#         if (verbose) cat(paste('iter=', niter), '\t',
#                          paste('param.error=', epsilon), '\t',
#                          paste('loglik=', LLnew), '\n');  
#         
#         #print(niter)
#         #print(zeta)
#         #print(SigmaPhi)
#         #print(SigmaE)
#         #print(LLnew)
#     }
#     
#     #Variance-covariance matrix
#     hessian <- hessian(loglik, theta, n2, Y, X, ID, n1s, 3)
#     var_cov <- solve(-hessian)
#     
#     param <- list(theta = list(zeta = zeta, 
#                                SigmaE = SigmaE, 
#                                SigmaPhi = SigmaPhi),
#                   loglik = LLnew, eps = epsilon, iter = niter,
#                   var_cov = var_cov)
#     return(param) 
#     
# }
