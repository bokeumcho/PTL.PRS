#' @title Function to convert p-values to correlation via the t-statistic
#' 
#' @references
#' Mak, T. S. H., Porsch, R. M., Choi, S. W., Zhou, X., & Sham, P. C. (2017).
#' Polygenic scores via penalized regression on summary statistics.
#' Genetic Epidemiology, 41(6), 469–480. \doi{10.1002/gepi.22050}
#'
#' See also the original lassosum package: \url{https://github.com/tshmak/lassosum}

#' @param n Sample size
#' @param p Vector of p-values
#' @param sign A vector giving the sign of the correlations (e.g. the log odds ratios)
#' @param min.n The minimum sample size to be considered a valid p-value
#' @return A vector of correlations
p2cor <- function(p, n, sign=rep(1, length(p)), min.n=max(n, 30)/10) {
  
  stopifnot(length(n)==1 || length(n) == length(p))
  stopifnot(length(p) == length(sign))
  
  t <- sign(sign) * qt(p/2, df=n-2, lower.tail=F)
  invalid <- n < min.n
  if(any(invalid)) {
    warning(paste(sum(invalid), "statistics has n < ", min.n, "and are coded as NA."))
  }
  t[invalid] <- NA
  
  return(t / sqrt(n - 2 + t^2))
  
}