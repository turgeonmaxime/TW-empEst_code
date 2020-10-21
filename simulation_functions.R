#' Colour palette for plotting results
#' @param n number of colours in the palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#' Fit empirical estimator to largest roots of permuted data
#'
#' This function uses the method of moments to fit the location-scale
#' Tracy-Widom family to the permuted values
#' 
#' @param null_dist vector of largest roots
fit_heuristic <- function(null_dist) {
  # Use method of moments
  muTW <- -1.2065335745820
  sigmaTW <- sqrt(1.607781034581)
  
  muS <- mean(log(null_dist))
  sigmaS <- stats::sd(log(null_dist))
  
  sigma1 <- sigmaS/sigmaTW
  mu1 <- muS - sigma1 * muTW
  
  return(c(mu1, sigma1))
}

# These functions extend the distribution functions from RMTstat
# to our location-scale family
#
#' @param x value at which to evaluate the density or distribution function
#' @param n number of observations
#' @param mu location parameter
#' @param sigma scale parameter
#' @param beta the order parameter (1, 2, or 4)
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].
dtw_ls <- function(x, mu, sigma, beta = 1, log = FALSE) {
  x1 <- (x - mu)/sigma
  # values for dtw are only available for the 
  # interval -10 to 6
  x1 <- as.numeric(x1 >= -10)*x1 - 10*as.numeric(x1 < -10)
  x1 <- as.numeric(x1 <= 6)*x1 + 6*as.numeric(x1 > 6)
  return(RMTstat::dtw(x1, beta, log)/sigma)
}

ptw_ls <- function(x, mu, sigma, beta = 1, lower.tail = TRUE, log = FALSE) {
  x1 <- (x - mu)/sigma
  # values for dtw are only available for the 
  # interval -10 to 6
  x1 <- as.numeric(x1 >= -10)*x1 - 10*as.numeric(x1 < -10)
  x1 <- as.numeric(x1 <= 6)*x1 + 6*as.numeric(x1 > 6)
  return(RMTstat::ptw(x1, beta, lower.tail, log.p = log))
}

rtw_ls <- function(n, mu, sigma, beta=1) {
  x <- RMTstat::rtw(n, beta)
  x1 <- x*sigma + mu
  return(x1)
}

# Generate a Wishart matrix, even in the singular setting
#
#' @param n number of observations
#' @param df degrees of freedom, which corresponds to the dimension of the corresponding normal random vector
#' @param Sigma scale matrix
#' @param shrink logical; if TRUE, use linear shrinkage estimate of the covariance instead of sample covariance
rWishart <- function(n, df, Sigma, shrink = FALSE) {
  list_mvrnorm <- replicate(n, mvtnorm::rmvnorm(df, sigma = Sigma,
                                                method = "chol"),
                            simplify = FALSE)
  list_Wishart <- if (shrink) {
    lapply(list_mvrnorm, function(X) shrink_est(X))
  } else lapply(list_mvrnorm, function(X) crossprod(X))
  
  return(list_Wishart)
}

# Solve determinantal equation 
#
#' @param A,B random matrices appearing in the determinantal equation
#' @param rankB rank of the matrix B
compute_largest_root <- function(A, B, rankB) {
  # rankB <- corpcor::rank.condition(B)$rank
  eigen_B <- mgcv::slanczos(B, k = rankB)
  eigen_B_values <- pmax(eigen_B$values, 0)
  nonzero_ind <- which(eigen_B_values > 0)
  transf <- eigen_B$vectors[,nonzero_ind] %*% diag(1/sqrt(eigen_B_values[nonzero_ind]))
  A_prime <- crossprod(transf, A %*% transf)
  
  rank_Ap <- corpcor::rank.condition(A_prime)$rank
  eigen_Ap <- mgcv::slanczos(A_prime, k = 1)
  largest_root <- max(eigen_Ap$values)
  
  return(sqrt(largest_root))
}

# Generate a Monte Carlo estimate of the largest root distribution
#
#' @param p,m,q parameters of the double Wishart problem
#' @param B number of Monte Carlo samples
#' @param rho parameter of the underlying exchangeable structure of the scale matrix Sigma
#' @param mc.cores number of cores to use for parallel computations
#' @param shrink logical; if TRUE, use linear shrinkage estimate of the covariance instead of sample covariance
generate_largestRootDist <- function(p, B, rho = 0, m = 4, q = 96, 
                                     mc.cores = 1, shrink = FALSE) {
  Sigma <- diag(nrow = p, ncol = p)
  if (rho != 0) {
    Sigma[upper.tri(Sigma)] <- Sigma[lower.tri(Sigma)] <- rho
  }
  
  Adist <- rWishart(B, df = m, Sigma = Sigma, shrink = FALSE)

  dist <- parallel::mclapply(Adist, function(A) {
    res <- mvtnorm::rmvnorm(q, sigma = Sigma, method = "chol")
    if (shrink) {
      B <- shrink_est(res)
      largestRoot <- compute_largest_root(A, B, rankB = p)
    } else {
      svdRes <- corpcor::fast.svd(res)
      rankVr <- corpcor::rank.condition(res)$rank
      eigVecVr <- svdRes$v[, 1:rankVr]
      eigValVrInv <- 1/svdRes$d[1:rankVr]
      Xp <- eigVecVr %*% diag(eigValVrInv)
      C <- crossprod(Xp, A %*% Xp)
      svdC <- corpcor::fast.svd(C)
      Xpp <- svdC$u
      singWeights <- Xp %*% Xpp
      largestRoot <- max(crossprod(singWeights,  A %*% singWeights)) 
    }
  }, mc.cores = mc.cores) 
  
  dist <- simplify2array(dist)

  return(dist)
  
}

# Linear shrinkage covariance estimate
#
#' @param res matrix of data
shrink_est <- function(res) {
  # port of matlab code from http://www.econ.uzh.ch/faculty/wolf/publications.html#9
  # Ledoit, O. and Wolf, M. (2004).
  # Honey, I shrunk the sample covariance matrix.
  # Journal of Portfolio Management 30, Volume 4, 110-119.
  p <- ncol(res); n <- nrow(res)
  
  # Compute sample covariance matrix using the de-meaned returns
  sample <- crossprod(scale(res, center = TRUE,
                            scale = FALSE))/n
  
  # Compute prior
  var <- matrix(diag(sample), ncol = 1)
  sqrtvar <- sqrt(var)
  tmpMat <- matrix(rep(sqrtvar, p), nrow = p)
  rBar <- (sum(sum(sample / (tmpMat * t(tmpMat)))) - p) / (p * (p - 1))
  prior <- rBar * tmpMat * t(tmpMat)
  diag(prior) <- var
  
  # What is called pi-hat
  y <- res^2
  phiMat <- crossprod(y) / n - 2 * crossprod(res) * sample / n + sample^2
  phi <- sum(phiMat)
  
  # What is called rho-hat
  term1 <- crossprod(res^3, res) / n
  help <- crossprod(res)/n
  helpDiag <- matrix(diag(help), ncol = 1)
  term2 <- matrix(rep(helpDiag, p), ncol = p, byrow = FALSE) * sample
  term3 <- help * matrix(rep(var, p), ncol = p, byrow = FALSE)
  term4 <- matrix(rep(var, p), ncol = p, byrow = FALSE) * sample
  thetaMat <- term1 - term2 - term3 + term4
  diag(thetaMat) <- 0
  rho <- sum(diag(phiMat)) + rBar * sum(sum(tcrossprod(1 / sqrtvar, sqrtvar) * thetaMat))
  
  # What is called gamma-hat
  gamma <- norm(sample - prior, "F")^2
  
  # Compute shrinkage constant
  kappa <- (phi - rho) / gamma
  shrinkage <- max(0, min(1, kappa / n))
  
  # Compute the estimator
  sigma <- shrinkage * prior + (1 - shrinkage) * sample
  sigma <- n*sigma
  
  return(sigma)
}
