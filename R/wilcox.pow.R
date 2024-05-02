#' Compute Power for the Wilcoxon and T Tests
#'
#' Computes an approximation to the power of the Wilcoxon signed rank or
#' two-sample Wilcoxon rank sum test and one- or two-sample t-tests
#'
#' @details
#' Uses an approximations given in Lehman, Nonparametrics: Statistical Methods
#' Based on Ranks, for the power of the Wilcoxon signed rank and two-sample
#' Wilcoxon test.  Assumes the data on a transformed scale is normally
#' distributed with constant variances.  The difference \code{del} is a
#' location shift on this transformed scale, expressed as a fraction of the
#' common standard deviation.
#'
#' The critical value from the t-distribution with n1-1 (one-sample) or n1+n2-2
#' (two-sample) degrees of freedom can be used instead of the standard normal
#' by specifying \code{tdist=TRUE}.  However, the rejection probability is
#' still computed using the normal distribution.
#'
#' @param n1 sample size in group 1
#' @param n2 Sample size in group 2. If NULL, assumes a 1-sample test.
#' @param del difference as a fraction of the population standard deviation
#' @param alpha2 two-sided significance level
#' @param tdist if TRUE, uses t distribution for the critical value; if FALSE,
#' uses a normal critical value
#'
#' @return
#' A vector of length two giving the power for the Wilcoxon (first component)
#' and t (second component) tests.
#'
#' @keywords design
#'
#' @examples
#' wilcox.pow(24, 24, 1)
#' wilcox.pow(24, 10, 1)
#' wilcox.pow(10, NULL, 0.5) ## one-sample test
#'
#' @export

wilcox.pow <- function(n1, n2 = NULL, del, alpha2 = 0.05, tdist = FALSE) {
  # del=difference in # standard deviations
  # two-sided test
  # del/2 factor is because variance of difference is 2sigma^2 and
  # normal density has a factor of 2pi
  df <- if (is.null(n2))
    n1 - 1 else n1 + n2 - 2
  q1 <- if (tdist)
    qt(1 - alpha2 / 2, df = df) else qnorm(1 - alpha2 / 2)
  if (is.null(n2)) {
    c(wilcox = pnorm((n1 * (n1 - 1) / 2 + n1 / sqrt(2)) * del / sqrt(pi) /
                       sqrt(n1 * (n1 + 1) * (2 * n1 + 1) / 24) - q1),
      t = pnorm(sqrt(n1) * del - q1))
  } else {
    c(wilcox = pnorm(sqrt(12 * n1 * n2 / ((n1 + n2 + 1) * pi)) * del / 2 - q1),
      t = pnorm(sqrt(n1 * n2 / (n1 + n2)) * del - q1))
  }
}
