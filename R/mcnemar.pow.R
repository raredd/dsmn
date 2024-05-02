#' Calculate Power for McNemar's Test
#'
#' Given a specified alternative distribution, this function computes the exact
#' power for NcNemar's Test
#'
#' @details
#' When two binary (for convenience, called positive vs. negative) factors are
#' measured on the same set of subjects, McNemar's test can be used to test
#' whether the marginal proportions positive for the two factors are the same.
#' For example, this could be of interest if the same quantity is measured
#' before and after treatment, or if two different assays supposedly measuring
#' the same quantity are being compared.
#'
#' For the values of the two factors, let p11=the proportion of the \code{n}
#' subjects positive for both, p22= the proportion negative for both, p12=the
#' proportion positive for factor 1 and negative for 2, and p21= the proportion
#' negative for 1 and positive for 2, where then p11+p22+p21+p12=1.  In terms
#' of these parameters, the input quantities are \code{pdisc} = p12+p21 and
#' \code{delta} = p12-p21.  If an alternative with p21>p12 is of interest, then
#' \code{delta} should be negative.
#'
#' McNemar's test is a conditional test of H0: p12=p21 (or equivalently the
#' marginal probabilities p11+p12 = p11+p21) based on the discordant pairs.
#' For the exact conditional test of one-sided size \code{alpha}, this function
#' uses the marginal distribution of the number of discordant pairs under the
#' specified distribution and the conditional distribution of the number
#' positive for factor 1 and negative for factor 2 under the null and the
#' alternative to determine the critical region of the test (under the null)
#' and the exact unconditional power (under the specified alternative).
#'
#'
#' @param n The total sample size
#' @param pdisc The combined proportion of discordant pairs
#' @param delta The difference in the discordant cell rates (factor 1 pos and
#' factor 2 neg minus factor 1 neg and factor 2 pos.
#' @param alpha The one-sided type I error of the test
#' @return A vector giving the values of \code{n}, p12, p21, the attained exact
#' size, the exact power, and Mietinen's approximation (Biometrics, 1968, p343)
#' to the power.
#'
#' @author Bob Gray
#'
#' @keywords design htest
#'
#' @examples
#' mcnemar.pow(100, 0.9, 0.2, delta = 0.095)
#'
#' @export mcnemar.pow

mcnemar.pow <- function(n, pdisc, delta, alpha = 0.025) {
  # pdisc = Combined probability of discordant pairs
  # delta = difference in success probabilities (ie between factor 1 and 2)
  # if the 4 cells in the 2x2 table have multinomial probabilities p11, p12,
  # p21, p22 then pdisc=p12+p21 and delta=p12-p21
  p12 <- (pdisc + delta) / 2
  p21 <- pdisc - p12
  if (min(p12, p21) < 0)
    stop('Incompatible probabilities specified')
  x <- 1:n
  xp <- dbinom(x, n, p12 + p21)
  size <- 0
  pow <- 0
  p <- max(p12, p21) / (p12 + p21)
  for (i in x) {
    z <- qbinom(1 - alpha, i, 0.5)
    if (z < i) {
      size <- size + (1 - pbinom(z, i, 0.5)) * xp[i]
      pow <- pow + (1 - pbinom(z, i, p)) * xp[i]
    }
  }
  # formula zb based on straightforward asymptotics
  # formula zb3 given by Mietinen as the straightforward version
  # formula zb3 from Mietinen Bcs, 1968, p.343 - given as more accurate
  zb2 <- (qnorm(1 - alpha) - abs(delta) / sqrt(pdisc / n)) /
    sqrt(1 - delta ^ 2 * (3 + pdisc) / (4 * pdisc ^ 2))
  c(n = n, p12 = p12, p21 = p21, size = size, power = pow, approx = 1 - pnorm(zb2))
}
