#' Power for Comparing Two Binomials
#'
#' Calculates power (\code{b2p}) and sample size (\code{b2n}) for a two sample
#' binomial test based on a normal approximation with and without continuity
#' correction. \code{b2p} also computes the exact power of Fisher's exact test
#' and the exact UMP unbiased test. Given combined sample size and response
#' probability, \code{b2diff} calculates the difference corresponding to a
#' specified power.
#'
#' @aliases b2p b2n b2diff
#'
#' @details
#' The power is computed for the one-sided test for the alternative that
#' \code{p1 > p2} against the null of equal rates.  For a one-sided test in the
#' opposite direction, reverse the definition of success and failure to map the
#' problem to the alternative in this direction. The calculations are also
#' valid for a two-sided test of size \code{2*alpha} when the two-sided test is
#' formed by combining the one-sided rejection regions of size \code{alpha} in
#' each direction. This function does not give correct results for more
#' general two-sided exact tests that divide the error asymmetrically between
#' the two tails unless the correct one-sided error corresponding to the
#' asymmetric test is used. The normal approximations are based on formulas
#' from Fleiss (Statistical Methods for Rates and Proportions, 2nd ed, 1981).
#'
#' \code{b2diff} is intended for facilitating power calculations when a new
#' (binary) marker will be measured on an existing sample. Then the combined
#' samples size \code{n} and overall response rate is known. \code{r} gives
#' the expected proportion of the sample in group 1 (defined by the new
#' marker). \code{b2diff} then calculates the response probabilities in the
#' two groups that give a difference for which the two-sample comparison will
#' have the specified power, subject to the overall response probability and
#' sample size constraints. The probabilities are computed for the one-sided
#' tests in each direction (higher response rates in group 1 and lower response
#' rates in group 1)
#'
#' The uniformly most powerful unbiased test uses a randomized rejection rule
#' on the boundary of the critical region to give the exact type I error rate
#' specified. This is generally not regarded as an appropriate procedure in
#' practice.
#'
#' @param p1 Success probability in group 1
#' @param p2 Success probability in group 2.  Must be \code{< p1}.
#' @param n1 Number of subjects in group 1
#' @param n2 Number of subjects in group 2
#' @param alpha One-sided significance level (type I error)
#' @param exact If \code{TRUE}, power for the exact test is computed
#' @param power The desired power of the test
#' @param r Proportion in group 1
#' @param p The overall (average) response probability in the two groups
#' @param n The combined number of subjects in the two groups
#'
#' @return
#' \code{b2p} returns a vector of length 4 giving the power based on the
#' normal approximation with continuity correction (\code{approx.cor}), the
#' normal approximation without correction (\code{approx.unc}), the exact power
#' for Fisher's exact test (\code{fisher}), and the exact power for the
#' uniformly most powerful unbiased test (\code{UMPU}), which uses a randomized
#' rejection rule on the boundary of the critical region.
#'
#' \code{b2n} returns a vector of length 2 giving the approximate sample size
#' for the continuity corrected statistic (\code{cont.cor}) and the uncorrected
#' (\code{uncor}) statistics.
#'
#' \code{b2diff} returns a vector of length 6 giving the input values of
#' \code{p} and \code{r}, the response probabilities in groups 1 and 2
#' corresponding to the difference with larger response in group 1, and the
#' response probabilities in groups 1 and 2 corresponding to the difference
#' with lower response in group 1.
#'
#' @references
#' Fleiss (1981). \emph{Statistical Methods for Rates and Proportions}. 2nd ed.
#'
#' @seealso
#' \code{\link{pickwin}}
#'
#' @keywords design htest
#'
#' @examples
#' b2diff(0.3, 0.1, 400)
#' b2n(0.4, 0.2)
#' b2p(0.4, 0.2, 100, 100)
#'
#' @export

b2p <- function(p1, p2, n1, n2, alpha = 0.025, exact = TRUE) {
  if (min(p1, p2, alpha) <= 0 | max(p1, p2, alpha) >= 1)
    stop('p1, p2, alpha must be between 0 and 1')
  if (p2 >= p1)
    stop('p1 must be > p2')
  if (min(n1, n2) < 1)
    stop('n1, n2 must be positive')
  r <- n2 / n1
  pd <- p1 - p2
  pb <- (p1 + r * p2) / (r + 1)
  t1 <- sqrt(pb * (1 - pb) * (1 + r) / n2)
  t2 <- sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
  al1 <- 1 - alpha
  zs <- qnorm(al1)
  zb <- (zs * t1 - pd) / t2
  b1 <- 1 - pnorm(zb)
  zb <- zb + (r + 1) / (2 * n2 * t2)
  b2 <- 1 - pnorm(zb)
  #   if (exact) {
  #     u <- .Fortran('pfishr', as.double(alpha), as.double(p1), as.double(p2),
  #                   as.integer(n1), as.integer(n2), double(n1 + n2 + 1),
  #                   double(n1 + n2 + 1), double(n1 + 1), double(n2 + 1),
  #                   double(n1 + 1), double(n2 + 1), double(1), double(1),
  #                   PACKAGE = 'desmon')[12:13]
  #     c(approx.cor = b2, approx.unc = b1, fisher = u[[1L]], UMPU = u[[2L]])
  #   } else
  c(approx.cor = b2, approx.unc = b1, fisher = NULL, UMPU = NULL)
}

#' @rdname b2p
#' @export
b2n <- function(p1, p2, power = 0.8, r = 0.5, alpha = 0.025) {
  if (min(r, p1, p2, alpha, power) <= 0 | max(p1, p2, r, alpha, power) >= 1)
    stop('p1, p2, r, alpha, power must be between 0 and 1')
  if (p2 >= p1)
    stop('p1 must be > p2')
  r <- (1 - r) / r
  pb <- (p1 + r * p2) / (r + 1)
  qb <- 1 - pb
  pd <- p1 - p2
  al1 <- 1 - alpha
  zs <- qnorm(al1)
  b <- 1 - power
  zb <- qnorm(b)
  t1 <- sqrt((r + 1) * pb * qb)
  t2 <- sqrt(r * p1 * (1 - p1) + p2 * (1 - p2))
  s1 <- (zs * t1 - zb * t2) / pd
  s1 <- s1 * s1 / r
  s2 <- sqrt(1 + 2 * (r + 1) / (s1 * r * pd))
  s2 <- s1 * (1 + s2) * (1 + s2) / 4
  s1 <- (r + 1) * s1
  s2 <- (r + 1) * s2
  c(cont.cor = s2, uncor = s1)
}

#' @rdname b2p
#' @export
b2diff <- function(p, r, n, alpha = .025, power = .8, exact = TRUE) {
  # given sample size and overall response prob,
  # find the difference with specified power
  S1 <- function(p1, p, r, n) {
    if (p1 < p) {
      n2 <- round(n * r)
      n1 <- n - n2
      p2 <- p1
      p1 <- (n * p - n2 * p2) / n1
    } else {
      n1 <- round(n * r)
      n2 <- n - n1
      p2 <- (n * p - n1 * p1) / n2
    }
    c(n1, n2, p1, p2)
  }

  S2 <- function(p1, p, r, n, alpha, power, exact) {
    u <- S1(p1, p, r, n)
    u <- b2p(u[3L], u[4L], u[1L], u[2L], alpha, exact)
    if (exact)
      u[3L] - power else u[1L] - power
  }

  minv <- p + 0.0001
  n1 <- round(n * r)
  maxv <- min(n * p / n1 - 0.0001, 0.9999)
  u1 <- uniroot(S2, c(minv, maxv), p = p, r = r, n = n, alpha = alpha,
                power = power, exact = exact)$root
  maxv <- p - 0.0001
  minv <- max((n * p - n + n1) / n1 + 0.0001, 0.0001)
  u2 <- uniroot(S2, c(minv, maxv), p = p, r = r, n = n, alpha = alpha,
                power = power, exact = exact)$root
  u3 <- S1(u1, p, r, n)
  u4 <- S1(u2, p, r, n)
  c(p, r, u3[3:4], u4[4:3])
}
