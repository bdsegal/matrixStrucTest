# Testing using Pearson and Filon's results -- Steiger's method

sigmaRhoFun <- function(j, k, h, m, A) {
  #' Sub-routine to compute OLS estimates of covariance between correlations
  #'
  #' This sub-routine is called by \code{X2fun}.
  #' @param j First index
  #' @param k Second index
  #' @param h Third index
  #' @param m Fourth index
  #' @param A Correlation matrix

  0.5 * (
  (A[j,h] - A[j,k]*A[k,h]) * (A[k,m] - A[k,h]*A[h,m]) +
  (A[j,m] - A[j,h]*A[h,m]) * (A[k,h] - A[k,j]*A[j,h]) +
  (A[j,h] - A[j,m]*A[m,h]) * (A[k,m] - A[k,j]*A[j,m]) +
  (A[j,m] - A[j,k]*A[k,m]) * (A[k,h] - A[k,m]*A[m,h]))
}

sigmaZFun <- function(s, t, index, A, Sigma) {
  #' Sub-routine to get covariance of z-transformed variables
  #'
  #' This sub-routine is called by \code{X2fun}.
  #' @param s First index
  #' @param t Second index
  #' @param index Matrix with two columns with index pairs given by rows
  #' @param A Correlation matrix
  #' @param Sigma Variance-covariances of correlation matrix A
  Sigma[s, t] / ((1 - A[index[s, 1], index[s, 2]]^2) * 
                 (1 - A[index[t, 1], index[t, 2]]^2))
}


X2Fun <- function(data, group_list, corMethod = "spearman"){
  #' Goodness-of-fit chi-squared statistic described by Steiger (1980).
  #'
  #' This function computed the goodness-of-fit chi-squared statistic described by Steiger (1980).
  #' Reference:
  #' Steiger, J. H. (1980). Tests for comparing elements of a correlation matrix. Psychological Bulletin, 87(2) 245-251.
  #' @param data Data frame
  #' @param group_list List of column indices of A for each group
  #' @param corMethod Type of correlations; passed to cor()
  #' @importFrom stats complete.cases cor pchisq t.test
  #' @export

  data <- data[complete.cases(data), ]
  N <- nrow(data)
  
  A <- cor(data, method = corMethod)

  p <- ncol(A)
  K <- length(group_list)
  
  if(sum(sapply(group_list,length)) > p){
    warning("Some variables may be assigned to multiple clusters")
  } else if(sum(sapply(group_list,length)) < p){
    warning("Some variables not assigned to a cluster")
  }
  
  # create delta matrix
  Delta <- matrix(nrow=p,ncol=p)
  for (i in 1:p){
    for (j in 1:p){
      Delta[i,j] <- deltaSub(i, j, group_list)
    }
  }
  
  # vectors of first and second indices in correlation matrix
  first <- NULL
  second <- NULL
  for (j in 2:p) {
  for (i in 1:(j-1)) {
    first <- c(first, i)
    second <- c(second, j)
    }
  }
  
  index <- cbind(first, second)
  nInd <- nrow(index)

  L <- cbind(1, Delta[upper.tri(Delta)])
  rho <- A[upper.tri(A)]

  # get initial OLS estimate and put in matrix for easy processing
  rhoHatOLSTri <- matrix(NA, nrow = nrow(A), ncol = ncol(A))
  rhoHatOLSTri[upper.tri(rhoHatOLSTri)] <- L %*% solve(crossprod(L, L), crossprod(L, rho))
  rhoHatOLSTri[lower.tri(rhoHatOLSTri)] <- t(rhoHatOLSTri)[lower.tri(rhoHatOLSTri)]
  diag(rhoHatOLSTri) <- 1

  # get SigmaOLS
  SigmaOLS <- matrix(0, nrow = nInd, ncol = nInd)
  for(s in 1:(nInd -1)) {
    for (t in (s):nInd) {
      SigmaOLS[s,t] <- sigmaRhoFun(index[s,1], index[s,2], 
                                   index[t,1], index[t,2], rhoHatOLSTri)
    }
  }
  SigmaOLS[lower.tri(SigmaOLS)] <- t(SigmaOLS)[lower.tri(SigmaOLS)]

  rhoHatGLS <- L %*% solve(crossprod(L, solve(SigmaOLS, L)), crossprod(L, 
               solve(SigmaOLS, rho)))

  # get covariance for Z transformed variables
  SigmaZ <- matrix(0, nrow = nInd, ncol = nInd)
  for (s in 1:(nInd - 1)) {
    for (t in s:nInd) {
      SigmaZ[s, t] <- sigmaZFun(s, t, index, A, SigmaOLS)
    }
  }
  SigmaZ[lower.tri(SigmaZ)] <- t(SigmaZ)[lower.tri(SigmaZ)]

  X2 <- (N - 3) *t(atanh(rho) -atanh(rhoHatGLS)) %*% 
        solve(SigmaZ, atanh(rho) - atanh(rhoHatGLS))
  df <- p*(p-1)/2 - 2

  result <- c(X2, df, pchisq(X2, df = df, lower.tail = FALSE))
  names(result) <- c("X2", "df", "pval")

  return(result)
}
