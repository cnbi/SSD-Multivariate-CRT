################### Loglikelihood function EM Multilevel ####################

loglik <- function(theta, n2, Y, X, ID, n1s, outcomes){
    
    if (outcomes == 2) {
        beta1 <- theta[1:2]
        beta2 <- theta[3:4]
        sphi11 <- theta[5]
        sphi12 <- theta[6]
        sphi22 <- theta[7]
        se11 <- theta[8]
        se12 <- theta[9]
        se22 <- theta[10]
        SigmaPhi <- matrix(c(sphi11, sphi12, sphi12, sphi22), 2, 2)
        SigmaE <- matrix(c(se11, se12, se12, se22), 2, 2)
        InvS2Phi <- solve(SigmaPhi)
        InvS2E <- solve(SigmaE)
        
        temp <- 0
        for (j in 1:n2) {
            Yj <- Y[ID == j, , drop = FALSE]
            Xj <- X[ID == j, , drop = FALSE]
            residj <- Yj - cbind(Xj %*% beta1, Xj %*% beta2)
            obs <- c(t(residj))
            tm1 <- (n1s[j] - 1) * log(det(SigmaE)) + log(det(SigmaE + n1s[j] * SigmaPhi))
            InvSS2 <- solve(SigmaE + n1s[j] * SigmaPhi) - InvS2E
            Invj <- kronecker(diag(nrow = n1s[j]),InvS2E) + kronecker(matrix(1, n1s[j], n1s[j]), InvSS2)/n1s[j]
            tm2 <- c(t(obs) %*% Invj %*% obs)
            temp <- temp - (tm1 + tm2)/2
        }
        
    } else if (outcomes == 3) {
        beta1 <- theta[1:2]
        beta2 <- theta[3:4]
        beta3 <- theta[5:6]
        
        sphi11 <- theta[7]
        sphi21 <- theta[8]
        sphi22 <- theta[9]
        sphi31 <- theta[10]
        sphi32 <- theta[11]
        sphi33 <- theta[12]
        
        se11 <- theta[13]
        se21 <- theta[14]
        se22 <- theta[15]
        se31 <- theta[16]
        se32 <- theta[17]
        se33 <- theta[18]
        
        SigmaPhi <- matrix(c(sphi11, sphi21, sphi31, sphi21, sphi22, sphi32, sphi31, sphi32, sphi33), 3, 3)
        SigmaE <- matrix(c(se11, se21, se31, se21, se22, se32, se31, se32, se33), 3, 3)
        InvS2Phi <- solve(SigmaPhi)
        InvS2E <- solve(SigmaE)
        
        temp <- 0
        for (j in 1:n2) {
            Yj <- Y[ID == j, , drop = FALSE]
            Xj <- X[ID == j, , drop = FALSE]
            residj <- Yj - cbind(Xj %*% beta1, Xj %*% beta2, Xj %*% beta3)
            obs <- c(t(residj))
            tm1 <- (n1s[j] - 1) * log(det(SigmaE)) + log(det(SigmaE + n1s[j] * SigmaPhi))
            InvSS2 <- solve(SigmaE + n1s[j] * SigmaPhi) - InvS2E
            Invj <- kronecker(diag(nrow = n1s[j]), InvS2E) + 
                kronecker(matrix(1, n1s[j], n1s[j]), InvSS2)/n1s[j]
            tm2 <- c(t(obs) %*% Invj %*% obs)
            temp <- temp - (tm1 + tm2)/2
        }
    }

    return(temp)
}
