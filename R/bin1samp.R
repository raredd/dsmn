#' Sample Size for One-sample Exact Binomial Tests
#'
#' Determines the minimum sample size for a one-sided, one-sample exact
#' binomial test with specified alpha (type I) and beta (type II) errors
#'
#' @details
#' Loops over sample sizes starting with \code{n.min}, determines the critical
#' value for the exact test of size alpha, and calculates the type II error.
#' If the type II error is larger than \code{beta}, increments the sample size
#' by 1 and tries again.
#'
#' \code{\link{print.bin1samp}} is a print method for the output of
#' \code{bin1samp}.
#'
#' @param p0 Null hypothesis response probability
#' @param pa Alternative hypothesis response probability
#' @param alpha Type I error rate
#' @param beta Type II error rate
#' @param n.min Minimum sample size considered
#'
#' @return
#' \code{bin1samp} returns a vector giving the minimum sample size
#' (\code{n}), the critical value \code{r} (reject if outcome is more extreme
#' than \code{r}), the null and alternative response probabilities (\code{p0}
#' and \code{pa}), and the type I and type II errors (\code{size} and
#' \code{type2}).)
#'
#' No value is returned by \code{print.bin1samp}, but it prints the input
#' \code{p0} and \code{pa}, the minimum sample size, the critical value
#' \code{r} for the test that rejects if the number of responses is \code{> r}
#' if \code{pa>p0} or that rejects if the number of responses is \code{< r} if
#' \code{pa<p0}, and the actual type I and type II error rates.
#'
#' @seealso
#' \code{\link{print.bin1samp}}; \code{\link{pickwin}}; \code{\link{rp21}};
#' \code{\link{twostg}}; \code{\link{simon}}
#'
#' @keywords design
#'
#' @examples
#' bin1samp(0.9, 0.95, n.min = 100)
#' bin1samp(0.1, 0.05, n.min = 100)
#'
#' @export

bin1samp <- function(p0, pa, alpha = 0.1, beta = 0.1, n.min = 20L) {
  if (p0 == pa)
    stop('p0 and pa are identical')
  b <- 1
  x <- round(p0 * n.min)
  n <- n.min - 1
  if (pa > p0) {
    while (b > beta) {
      n <- n + 1
      # determine cutoff: reject if X>x
      l <- x:n
      s <- 1 - pbinom(l, n, p0)
      sub <- s <= alpha
      x <- l[sub][1L]
      size <- s[sub][1L]
      b <- pbinom(x, n, pa)
    }
  } else if (pa < p0) {
    while (b > beta) {
      n <- n + 1
      # determine cutoff: reject if X<x
      l <- x:0
      s <- pbinom(l, n, p0)
      sub <- s <= alpha
      x <- l[sub][1L]
      size <- s[sub][1L]
      b <- 1 - pbinom(x, n, pa)
      x <- x + 1
    }
  }

  structure(
    c(n = n, r = x, p0 = p0, pa = pa, size = size, type2 = b),
    class = 'bin1samp'
  )
}

#' Print a summary of the output from bin1samp
#'
#' Prints a summary of the output from \code{\link{bin1samp}}.
#'
#' @details
#' Prints operating characteristics and stopping rules of
#' \code{\link{bin1samp}}.
#'
#' @param x An object (vector) of class \code{bin1samp}
#' @param ... Included for compatibility with the generic function
#'
#' @return
#' No value is returned - used for side-effect of printing to console.
#'
#' @seealso
#' \code{\link{bin1samp}}
#'
#' @export

print.bin1samp <- function(x, ...) {
  cat(paste0('             p0 = ', format(x[3L]), '\n'))
  cat(paste0('             pa = ', format(x[4L]), '\n'))
  cat(paste0('min sample size = ', x[1L],'\n\n'))
  if (x[3L] < x[4L])
    cat('reject H0 if # responses >', format(x[2L]), '\n\n')
  else cat('reject H0 if # responses <', format(x[2L]), '\n\n')
  cat(paste0('   type I error = ', format(x[5L]), '\n'))
  cat('  type II error =', format(x[6L]),'\n')
  invisible(NULL)
}
