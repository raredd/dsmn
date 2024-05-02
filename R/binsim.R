#' Simulation for Power of a Sequential Fisher's Exact Test
#'
#' Performs a simulation to estimate the power of Fisher's Exact under group
#' sequential monitoring.
#'
#' @details
#' Calls \code{\link{fisher.test}}. If the experimental group (2) has a higher
#' event rate than the control group (1), then 'less' is the appropriate
#' direction for a one-sided test.
#'
#' @param nsamp Number of simulated trials
#' @param n Total sample size for each trial
#' @param p1 Success probability for group 1 (control)
#' @param p2 Success probability for group 2
#' @param inf.times Information times of analyses (as fractions)
#' @param upper Upper boundary (p-value scale, so test rejects for values
#' smaller than the 'upper' boundary)
#' @param lower Lower boundary (p-value scale)
#' @param p.con Proportion randomized to control (group 1)
#' @param alt Direction of the alternative (one of greater, less or two.sided)
#'
#' @return A vector of length \code{nsamp}, giving the test result for each
#' simulated trial, coded 1 if the upper boundary was crossed, -1 if the lower
#' boundary was crossed, and coded 0 if neither was (only the first boundary
#' crossing is recorded).
#'
#' @author R Gray
#'
#' @seealso
#' \code{\link{b2p}}; \code{\link{fisher.test}}
#'
#' @keywords design
#'
#' @examples
#' binsim(20, 50, 0.1, 0.5, c(0.5, 1), c(0.05, 0.05), c(0.5, 0.5), alt = 'less')
#'
#' @export

binsim <- function(nsamp, n, p1, p2, inf.times, upper,
                   lower = rep(1.5, length(inf.times)), p.con = 0.5,
                   alt = 'less') {
  # upper and lower need to be in terms of one-sided p-values
  if (length(upper) != length(inf.times))
    stop('length upper != length inf.times')
  out <- rep(0, nsamp)
  for (i in 1:nsamp) {
    n1 <- c(0, round(p.con * inf.times * n))
    n2 <- c(0, round((1 - p.con) * inf.times * n))
    x1 <- cumsum(rbinom(length(inf.times), diff(n1), p1))
    x2 <- cumsum(rbinom(length(inf.times), diff(n1), p2))
    j <- 0
    while (j < length(inf.times) & out[i] == 0) {
      j <- j + 1
      u <- matrix(c(x1[j], x2[j], n1[j + 1L] - x1[j], n2[j + 1L] - x2[j]), 2L)
      u <- fisher.test(u, alternative = alt)$p.value
      if (u <= upper[j])
        out[i] <- 1
      if (u >= lower[j])
        out[i] <- -1
    }
  }
  out
}
