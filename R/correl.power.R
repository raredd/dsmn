#' Compute Power for Failure Time Comparisons in Correlative Studies
#'
#' Computes power or difference for the logrank test for a correlative study
#' comparing failure time distributions between groups defined by two levels of
#' a correlative marker when there is a known amount of information (number of
#' events)
#'
#' @details
#' Assumes exponential failures and a clinical trial setting with uniform
#' accrual over \code{acc.per} units of time plus an additional \code{add.fu}
#' time units of follow-up, so censoring is uniform on (\code{add.fu},
#' \code{add.fu + acc.per}).  Test statistic is the exponential Wald test using
#' asymptotic normality, which is asymptotically equivalent to the logrank
#' test.
#'
#' Only the variance under the alternative is used, consistent with
#' Rubinstein, Gail and Santner (1981) and \code{seqopr} and \code{srvpwr}.
#' \code{r1o2} can be either larger or smaller than 1.  The total failure
#' probability is fixed at \code{nd/n}, and hazard rates in the two groups are
#' computed to give the specified ratio subject to the total failure
#' probability constraint.  The \code{uniroot} function is called to solve for
#' the hazards, and \code{int} is the initial bracketing interval used in the
#' search.  Since the magnitude of the hazard functions are dependent on the
#' time units, it may sometimes be necessary to modify this interval.
#'
#' Given a specified ratio, \code{correl.power} computes the power.  Given a
#' specified power, \code{correl.ratio} computes the corresponding hazard ratio
#' (two ratios are computed, one for the difference in each direction).
#'
#' @aliases correl.power correl.ratio
#'
#' @param n Total number of cases
#' @param nd Total number of failures
#' @param r1o2 Ratio of the hazard in Group 1 over that in Group 2
#' @param p1 Proportion of sample in Group 1
#' @param acc.per Duration of accrual
#' @param add.fu Additional follow-up after accrual was completed
#' @param alpha2 Type I error rate (two-sided)
#' @param int Initial search interval for calculating hazard rates
#' @param power The desired power for the comparison
#' @param maxrat The maximum hazard ratio to consider in the search when
#' \code{greater=TRUE}. \code{1/maxrat} is the minimum ratio to consider when
#' \code{less=TRUE}
#' @param greater logical; if \code{TRUE}, find the ratio greater than 1 with
#' the specified \code{power}
#' @param less logical; if \code{TRUE}, find the ratio less than 1 with the
#' specified \code{power}
#'
#' @return
#' \code{correl.power} returns a vector giving the hazard rates in
#' groups 1 and 2 (\code{lambda1} and \code{lambda2}), the failure
#' probabilities in the two groups (\code{fp1} and \code{fp2}), the standard
#' error of the log of the ratio of the estimated exponential failure rates
#' (\code{se}), the total expected number of failures under the specified
#' failure rates (\code{n.fail}), the normal deviate corresponding to the type
#' II error (\code{zbeta}), and the power (\code{power}).
#'
#' \code{correl.ratio} returns a vector giving the total number of failures
#' (\code{n.fail}), the proportion of the sample expected to be in group 1
#' (\code{pr.g1}), and the hazard ratios (group 1/group 2) for which the test
#' will have the specified power in the directions with worse and better
#' outcome on group 1 (\code{r1o2.1worse} and \code{r1o2.1better}).
#'
#' @seealso \code{\link{powlgrnk}}
#'
#' @references
#' Rubinstein, Gail, Santner (1981). \emph{J Chron Dis} \strong{8}: 67-74.
#' Rubinstein, Gail, Santner (1981). \emph{J Chron Dis} \strong{34}: 469-479.
#'
#' @keywords design survival
#'
#' @examples
#' correl.power(355,172,1.54,.4,7,4)
#' correl.power(355,172,3.36,.92,7,4)
#' correl.power(355,172,1.67,.7,7,4)
#'
#' # check symmetry
#' correl.power(355,172,2/3,.5,7,4)
#' correl.power(355,172,1.5,.5,7,4)
#' correl.power(355,172,2/3,.25,7,4)
#' correl.power(355,172,1.5,.75,7,4)
#'
#' correl.ratio(1000,150,.2,2,3.635,int=c(.0001,20))
#'
#' @export correl.power

correl.power <- function(n, nd, r1o2, p1, acc.per, add.fu,
                         alpha2 = 0.05, int = c(0.001, 5)) {
  zc <- qnorm(1 - alpha2 / 2)
  q <- 1 - p1
  fpnull <- nd / n
  a <- add.fu
  b <- acc.per + add.fu # censoring uniform on (a,b)
  cpx <- function(l) (exp(-a * l) - exp(-b * l)) / (l * acc.per)
  fp <- function(l) {
    h <- r1o2 * l
    1 - q * cpx(l) - p1 * cpx(h) - fpnull
  }
  z <- uniroot(fp, int)
  h <- z[[1L]]
  l <- h * r1o2
  fp1 <- 1 - cpx(l)
  fp2 <- 1 - cpx(h)
  # use sa under null and alt, per RGS 1981 -- also, consistent with seqopr and srvpwr
  sa <- sqrt(1 / (n * p1 * fp1) + 1 / (n * q * fp2))
  zb <- (zc * sa - abs(log(r1o2))) / sa
  c(lambda1 = l, lambda2 = h, fp1 = fp1, fp2 = fp2, se = sa,
    n.fail = n * (fp1 * p1 + fp2 * q), zbeta = zb, power = 1 - pnorm(zb))
}

#' @export
correl.ratio <- function(n, nd, p1, acc.per, add.fu, power = 0.8, alpha2 = 0.05,
                         int = c(0.001, 5), maxrat = 20, greater = TRUE,
                         less = TRUE) {
  ## given power and information, find the ratio
  F1 <- function(rat, nd, p1, n, acc.per, add.fu, alpha2, power, intb) {
    correl.power(n, nd, rat, p1, acc.per, add.fu, alpha2, intb)[8L] - power
  }
  u1 <- if (greater) {
    uniroot(F1, c(1.01, maxrat), n = n, nd = nd, p1 = p1, acc.per = acc.per,
            add.fu = add.fu, alpha2 = alpha2, power = power, intb = int)$root
  } else NULL
  u2 <- if (less) {
    uniroot(F1, c(0.99, 1 / maxrat), n = n, nd = nd, p1 = p1, acc.per = acc.per,
            add.fu = add.fu, alpha2 = alpha2, power = power, intb = int)$root
  } else NULL

  c(n.fail = nd, pr.g1 = p1, r1o2.1worse = u1, r1o2.1better = u2)
}
