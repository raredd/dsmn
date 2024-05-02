#' Exact Binomial Confidence Limits
#'
#' Calculates exact confidence bounds on a single binomial probability
#'
#' @details
#' Calculates exact \code{1-(1-conf)/2} upper and lower confidence limits for a
#' binomial proportion.  The precision of the solution is the default in
#' \code{uniroot}.
#'
#' @param r number of successes
#' @param n total number of trials
#' @param conf confidence level
#'
#' @return
#' Returns a vector containing the two confidence limits. Returns 0 for
#' the lower limit if r=0, and returns 1 for the upper limit if r=n.
#'
#' @seealso
#' \code{\link{twocon}}
#'
#' @keywords htest
#'
#' @examples
#' binci(54, 88)
#'
#' @export binci

binci <- function(r, n, conf = 0.95) {
  if (r < 0 | r > n)
    stop('invalid value for r')
  if (conf <= 0 | conf >= 1)
    stop('must have 0 < conf <1')
  alpha <- (1 - conf) / 2
  ff <- function(p, r, n, alpha) 1 - pbinom(r - 1, n, p) - alpha
  pl <- if (r <= 0)
    0 else uniroot(ff, c(1e-8, 1 - 1e-8), r = r, n = n, alpha = alpha)$root
  ff2 <- function(p, r, n, alpha) pbinom(r, n, p) - alpha
  pu <- if (r >= n)
    1 else uniroot(ff2, c(1e-8, 1 - 1e-8), r = r, n = n, alpha = alpha)$root
  c(pl, pu)
}

#' One Sample Binomial Design
#'
#' Determines the critical value and power for an exact one sample binomial
#' test
#'
#' @details
#' Let \code{R} be the number of responses.  If \code{p0<pa}, determines the
#' smallest critical value \code{r} such that \code{P(R>=r|p0)<=alpha}, and
#' computes \code{P(R>=r|p)} for \code{p=p0,pa}.  If \code{p0>pa}, similar
#' calculations are done for rejecting in the opposite tail.
#'
#' @param n sample size
#' @param p0 The null response probability
#' @param pa The alternative response probability
#' @param alpha The one-sided type I error rate
#' @return A vector giving the critical value and the rejection probabilities
#' under the null and alternative.
#'
#' @seealso
#' \code{\link{bin1samp}}
#'
#' @keywords design
#'
#' @examples
#' onearm(30, 0.1, 0.3)
#'
#' @export onearm

onearm <- function(n, p0, pa, alpha = 0.1) {
  # single arm single stage evaluation
  u <- pbinom(0:n, n, p0) # P(<=)
  if (p0 > pa) {
    index <- max((1:(n + 1))[u <= alpha])
    c('rej if <=' = index - 1, size = u[index],
      power = pbinom(index - 1, n, pa))
  } else {
    u <- c(1, 1 - u[-n - 1]) # P(>=)
    index <- min((1:(n + 1))[u <= alpha])
    c('rej if >=' = index - 1, size = u[index],
      power = 1 - pbinom(index - 2, n, pa))
  }
}

#' Outcome Probabilities for Randomized Phase II Designs
#'
#' For a single stage randomized phase II study, computes the probability that
#' the difference in the number of responses is larger than a specified
#' critical value.
#'
#' @details
#' Computations assume two independent binomials R1 and R2, with R1 distributed
#' Binomal(n,p1) and R2 distributed Binomial(n,p2)
#'
#' @param n Planned number of subjects per arm
#' @param p1 Success probability for treatment 1
#' @param p2 Success probability for treatment 2
#' @param crit Maximum difference in number of responses that would be viewed
#' as an equivalent outcome
#'
#' @return
#' A vector giving P(R1>R2+crit), P(|R1-R2|<=crit) and P(R1<R2-crit).
#'
#' @seealso
#' \code{\link{pickwin}}
#'
#' @keywords design
#'
#' @examples
#' rp21(32, 0.2, 0.1)
#'
#' @export rp21

rp21 <- function(n, p1, p2, crit = 0) {
  # computes joint probs for two indep binom(n, p1), n(p2)
	probs <- outer(dbinom(0:n, n, p1), dbinom(0:n, n, p2))
	rw <- row(probs)
	cl <- col(probs)
	c('P(R1>R2+crit)' = sum(probs[rw > cl + crit]),
    'P(|R1-R2|<=crit)' = sum(probs[abs(rw - cl) <= crit]),
		'P(R1<R2-crit)' = sum(probs[rw < cl - crit]))
}

#' Two-Stage Randomized Phase II Designs
#'
#' Calculates operating characteristics of pick the winner rules for two stage
#' randomized phase II designs
#'
#' @details
#' Assumes that during the first stage of accrual \code{n1} patients are
#' randomized to each treatment (RX1 and RX2).  If more than \code{r1}
#' (\code{r3}) responses are seen on treatment RX1 (RX2), then an additional
#' \code{n2} patients will be accrued to treatment RX1 (RX2).  If more than
#' \code{r2} (\code{r4}) total responses are seen on RX1 (RX2), then RX1 (RX2)
#' will be declared active.
#'
#' Let X1 and X2 be the number of responses observed on treatments 1 and 2.
#' This function computes and prints P(X1<=r2)=P(RX1 declared inactive),
#' P(X2<=r4)=P(RX2 declared inactive), P(X2>r4 and X2>X1)=P(RX2 declared
#' superior), P(X1>r2 and X1>X2)=P(RX1 declared superior), and P(X1=X2 or
#' (X1<=r2 and X2<=r4))=P(neither declared superior).  Note that P(neither
#' declared superior) includes both the probability that both are declared
#' inactive and the probability that both are declared active and are equal.
#'
#' @param n1 first stage sample size (assumed to be the same for both
#' treatments)
#' @param n2 addition accrual during the second stage (assumed to be the same
#' for both treatments)
#' @param p1 True response probability on treatment 1
#' @param p2 True response probability on treatment 2
#' @param r1 The maximum number of responses to declare treatment 1 inactive in
#' the first stage
#' @param r2 The maximum number of total responses over both stages to declare
#' treatment 1 inactive
#' @param r3 The maximum number of responses to declare treatment 2 inactive in
#' the first stage
#' @param r4 The maximum number of total responses over both stages to declare
#' treatment 2 inactive
#'
#' @return
#' Returns a vector of length 9, whose components are (in order)
#' P(X1<=r1, X2>r4), P(r1<X1<=r2, X2>r4), P(X1>r2, X2>X1, X2>r4), P(X2<=r3,
#' X1>r2), P(r3<X2<=r4, X1>r2), P(X1>r2, X2>r4, X1>X2), P(X1<=r2, X2<=r4),
#' P(X1>r2, X2>r4, X1=X2), and the sum of the first 8 terms.
#'
#' @seealso
#' \code{\link{b2p}}
#'
#' @references
#' Simon R, Wittes RE and Ellenberg SS, 1985. Randomized phase II
#' clinical trials. \emph{Cancer Treat Rep} \strong{69} (12):1375--81.
#'
#' @keywords design
#'
#' @examples
#' pickwin(14, 18, 0.1, 0.2, 0, 3)
#'
#' @export pickwin

pickwin <- function(n1, n2, p1, p2, r1, r2, r3 = r1, r4 = r2) {
  # prob calculations for 2-stage pick the winner design
  # n1=# entered at first stage, n2=additional # entered at second stage
  # p1,p2 response probabilities on treatments 1 and 2
  # r1=max number responses stage 1 trt 1 and still stop
  # r2=max total number responses on trt 1 to be inactive
  # r3,r4 = corresponding values of r1,r2 for trt 2
  # prob( both are active and get the same # responses) included in p(neither
  #    better)
  if (n1 < 1 | n2 < 0 | r1 < 0 | r3 < 0 | r2 < 0 | r4 < 0 | p1 <= 0 | p2 <= 0 |
      r1 > n1 | r3 > n1 | r2 > n2 + n1 | r4 > n2 + n1 | p1 >= 1 | p2 >=1)
    stop('invalid arguments')
  x1 <- 0:n1
  x2 <- 0:n2
  w1 <- dbinom(x1, n1, p1)
  w2 <- dbinom(x1, n1, p2)
  w3 <- dbinom(x2, n2, p1)
  w4 <- dbinom(x2, n2, p2)
  # b1 = marginal dist of # responses in group 1
  u1 <- c(outer(x1[-(1:(r1 + 1))], x2, '+'))
  u2 <- c(outer(w1[-(1:(r1 + 1))], w3))
  b1 <- c(w1[1:(r1 + 1)], tapply(u2, u1, sum))
  # b2 = marginal dist of # responses in group 2
  u1 <- c(outer(x1[-(1:(r3 + 1))], x2, '+'))
  u2 <- c(outer(w2[-(1:(r3 + 1))], w4))
  b2 <- c(w2[1:(r3 + 1)], tapply(u2, u1, sum))
  # print(c(sum(b1), sum(b2)))

  cat('RX1 is declared inactive if <=', format(r1), 'responses are observed',
      ' in\nthe first', format(n1), 'cases or <=', format(r2), 'responses are',
      'observed in the first', format(n1 + n2), 'cases\n')
  cat('P(RX1 declared inactive) =', format(signif(sum(b1[1:(r2 + 1)]), 3)), '\n\n')
  cat('RX2 is declared inactive if <=', format(r3), 'responses are observed',
      ' in\nthe first', format(n1), 'cases or <=', format(r4), 'responses are',
      'observed in the first', format(n1 + n2), 'cases\n')
  cat('P(RX2 declared inactive) =', format(signif(sum(b2[1:(r4 + 1)]), 3)), '\n\n')
  a1 <- sum(b1[1:(r1 + 1)]) * sum(b2[(r4 + 2):length(b2)])
  a2 <- if (r2 > r1)
    sum(b1[(r1 + 2):(r2 + 1)]) * sum(b2[(r4 + 2):length(b2)]) else 0
  a3 <- 0
  for (i in (r2 + 2):(length(b1) - 1))
    a3 <- a3 + b1[i] * sum(b2[(max(i + 1, r4 + 2)):length(b2)])
  cat('prob RX2 declared better:', format(signif(a1 + a2 + a3, 3)), '\n')
  a4 <- sum(b2[1:(r3 + 1)]) * sum(b1[(r2 + 2):length(b1)])
  a5 <- if (r4 > r3)
    sum(b2[(r3 + 2):(r4 + 1)]) * sum(b1[(r2 + 2):length(b1)]) else 0
  a6 <- 0
  for (i in (r4 + 2):(length(b2) - 1))
    a6 <- a6 + b2[i] * sum(b1[(max(i + 1, r2 + 2)):length(b1)])
  cat('prob RX1 declared better:', format(signif(a4 + a5 + a6, 3)), '\n')
  a7 <- sum(b1[1:(r2 + 1)]) * sum(b2[1:(r4 + 1)])
  i <- (max(r2 + 2, r4 + 2)):length(b1)
  a8 <- sum(b1[i] * b2[i])
  cat('prob neither declared better:', format(signif(a7 + a8, 3)), '\n')
  c(a1, a2, a3, a4, a5, a6, a7, a8, a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8)
}

#' Confidence Interval on the Response Rate from a Two Stage Study
#'
#' Computes a confidence interval and several estimators of the response rate
#' using data from a phase II study with two-stage sampling
#'
#' @details
#' First \code{n1} patients are entered on the study.  If more than \code{r1}
#' responses are observed, then an additional \code{n2} patients are entered.
#' This function assumes that if the observed number of response \code{r < r1},
#' then only \code{n1} patients were entered.
#'
#' The estimators computed are the MLE (the observed proportion of responses),
#' a bias corrected MLE, and an unbiased estimator, which is sometimes
#' incorrectly described as the UMVUE.
#'
#' The confidence interval is based on the exact sampling distribution.
#' However, there is not a universally accepted ordering on the sample space in
#' two-stage designs.  The parameter \code{dp} can be used to modify the
#' ordering by weighting points within the sample space differently.
#' \code{dp=0} will give the Atkinson and Brown procedure, and \code{dp=1} will
#' order outcomes base on the MLE.  The Atkinson and Brown procedure orders
#' outcomes based solely on the number of responses, regardless of the number
#' cases sampled.  The MLE ordering defines as more extreme those outcomes with
#' a more extreme value of the MLE (the proportion of responses).  Other powers
#' of \code{dp}, such as \code{dp=1/2}, could also be used.  Let \code{R} be
#' the number of responses and \code{N=n1} if \code{R<=r1} and \code{N=n1+n2}
#' if \code{R>r1}.  In general, the outcomes that are more extreme in the high
#' response direction are those with \code{R/(N^dp) >= r/(n^dp)}, where
#' \code{r} and \code{n} are the observed values of \code{R} and \code{N}, and
#' the outcomes that are more extreme in the low response direction are those
#' with \code{R/(N^dp) <= r/(n^dp)}.
#'
#' @param n1 Number of cases entered during the first stage
#' @param n2 Number of additional cases to be entered during the second stage
#' @param r1 max number of responses that can be observed in the first stage
#' without continuing
#' @param r total number responses observed
#' @param conf two-sided confidence level (proportion) for the confidence
#' interval
#' @param dp Affects the ordering of outcomes within the sample space (see
#' below)
#'
#' @return
#' A vector with the lower confidence limit, the upper confidence
#' limit, the bias corrected MLE, the MLE, and the unbiased estimator.
#'
#' @seealso
#' \code{\link{binci}}
#'
#' @references
#' Atkinson and Brown (1985), BIOMETRICS 741-744.
#'
#' @keywords htest
#'
#' @examples
#' twocon(14, 18, 3, 4, dp = 0)
#' twocon(14, 18, 3, 4, dp = 1)
#'
#' @export twocon

twocon <- function(n1, n2, r1, r, conf = 0.95, dp = 1) {
  # confidence interval for two-stage phase II
  # n1 = # entered at first stage, n2=additional # entered at second stage
  # r1 = max number responses stage 1 trt 1 and still stop
  # r = total number responses observed
  if (n1 < 1 | n2 < 1 | r1 < 0 | r1 > n1 | r < 0 |
      r > n2 + n1 | conf <= 0 | conf >= 1)
    stop('invalid arguments')
  alpha <- (1 - conf) / 2
  x1 <- 0:n1
  x2 <- 0:n2
  u1 <- c(outer(x1[-(1:(r1 + 1))], x2, '+'))
  dbin2 <- function(p1, x1, x2, u1, n1, n2, r1) {
    w1 <- dbinom(x1, n1, p1)
    w3 <- dbinom(x2, n2, p1)
    u2 <- c(outer(w1[-(1:(r1 + 1))], w3))
    # ith component is P(R=i-1)
    c(w1[1:(r1 + 1)], tapply(u2, u1, sum))
  }
  n <- n1 + n2
  mle <- if (r > r1) r/n else r / n1
  # bias corrected estimator
  mm <- c((0:r1) / n1, (r1 + 1):n / n)
  ff <- function(p, x1, x2, u1, n1, n2, r1, mm, dbin2, mle)
    sum(dbin2(p, x1, x2, u1, n1, n2, r1) * mm) - mle
  pm <- if (r <= 0)
    0 else if (r >= n)
      1 else uniroot(ff, c(1e-8, 1 - 1e-8), x1 = x1, x2 = x2, u1 = u1, n1 = n1,
                     n2 = n2, r1 = r1, mm = mm, dbin2 = dbin2, mle = mle)$root
  # unbiased estimator
  ube <- if (r1 >= r) {
    r/n1
  } else {
    aa <- dhyper((r1 + 1):r, n1, n2, r)
    sum(((r1 + 1):r) * aa) / (n1 * sum(aa))
  }
  mm <- if (r > r1)
    c((n / n1) ^ dp * (0:r1),( r1 + 1):n) else c(0:r1,(n1 / n) ^ dp * ((r1 + 1):n))
  s1 <- mm >= r # for P(R>=r)
  s2 <- mm <= r
  ff2 <- function(p, x1, x2, u1, n1, n2, r1, s, dbin2, alpha)
    sum(dbin2(p, x1, x2, u1, n1, n2, r1)[s]) - alpha
  pl <- if (r <= 0)
    0 else uniroot(ff2, c(1e-8, 1 - 1e-8), x1 = x1, x2 = x2, u1 = u1, n1 = n1,
                   n2 = n2, r1 = r1, s = s1, dbin2 = dbin2, alpha = alpha)$root
  pu <- if (r >= n)
    1 else uniroot(ff2, c(1e-8, 1 - 1e-8), x1 = x1, x2 = x2, u1 = u1, n1 = n1,
                   n2 = n2, r1 = r1, s = s2, dbin2 = dbin2, alpha = alpha)$root
  c(lower = pl, upper = pu, bcmle = pm, mle = mle, unbiased = ube)
}
