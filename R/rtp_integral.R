
#' RTP p-value by integral approximation
#'
#' Computes RTP p-values using an integral approximation.
#'
#' @param p.values Numeric vector of p-values.
#' @param i Truncation point or vector of truncation points.
#'
#' @return Numeric RTP p-value or vector of RTP p-values.
#'
#' @export
#'
rtp_integral <- function(p.values, i) {
  v.stat.pval <- sapply(i, function(i) rtp::p.rtp(i, p.values))
  return(v.stat.pval)
}
