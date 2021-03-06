#' This function computes minimum square estimates for Fourier coefficients.
#'
#' @param f A time series to be smoothed.
#' @param m Order of the Fourier polynomial. Default value is computed using
#' the Sturge's rule.
#'
#' @return A vector with the fourier coefficients.
#' @import magrittr
#'
#' @examples
#' data(canada)
#'
#' coef_fourier(canada$ThePas_ts)
#' @export
coef_fourier <- function(f, m) {
  #m: ordem in the fourier approximation
  #f: time series to smooth

  if (missing(m)) m <-
      ceiling(1 + length(f) %>% log2())

  n <- length(f) #length of time series

  u <-  pi + 2 * pi * (((1:n) - 1) / n - 1)

  m_a <- matrix(0, nrow = 2 * m + 1, ncol = 2 * m + 1)
  b <- matrix(0, nrow = 2 * m + 1, ncol = 1)

  #filling A
  #a_0
  m_a[1, 1] <- n / 2
  for (k in 1:n) {
    m_a[1, 2:(m + 1)] <- m_a[1, 2:(m + 1)] + cos((1:m) * u[k])
    m_a[1, (m + 2):(2 * m + 1)] <- m_a[1, (m + 2):(2 * m + 1)] +
      sin((1:m) * u[k])
  }


  for (j in 1:m) {
    m_a[j + 1, 1] <- sum(cos(u) * j / 2)
    m_a[m + 1 + j, 1] <- sum(sin(u) * j / 2)
    for (k in 1:n) {
      m_a[j + 1, 2:(m + 1)] <- m_a[j + 1, 2:(m + 1)] +
        cos((1:m) * u[k]) * cos(j * u[k])
      m_a[j + 1, (m + 2):(2 * m + 1)] <- m_a[j + 1, (m + 2):(2 * m + 1)] +
        sin((1:m) * u[k]) *  cos(j * u[k])

      m_a[m + j + 1, 2:(m + 1)] <- m_a[m + j + 1, 2:(m + 1)] +
        cos((1:m) * u[k]) *  sin(j * u[k])
      m_a[m + j + 1, (m + 2):(2 * m + 1)] <-
        m_a[m + j + 1, (m + 2):(2 * m + 1)] + sin((1:m) * u[k]) * sin(j * u[k])
    }
  }

  #filling b
  b[1] <- sum(f)
  for (k in 1:n) {
    b[2:(m + 1)] <- b[2:(m + 1)] +
      f[k] * cos((1:m) * u[k])
    b[(m + 2):(2 * m + 1)] <- b[(m + 2):(2 * m + 1)] +
      f[k] * sin((1:m) * u[k])
  }

  #solving the linear system
  solve(m_a, b)
}
