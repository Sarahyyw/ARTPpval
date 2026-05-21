
#' ARTP p-value by ultra-fast interpolation
#'
#' Computes an ARTP p-value using a precomputed interpolation table.
#'
#' @param p.values Numeric vector of p-values.
#'
#' @return A list with estimated ARTP p-value, selected cutpoint,
#' observed ARTP statistic, and RTP p-values across truncation points.
#'
#' @export
artp_ufi <- function(p.values) {

  if (!is.numeric(p.values)) {
    stop("p.values must be numeric.")
  }

  if (any(is.na(p.values))) {
    stop("p.values cannot contain missing values.")
  }

  if (any(p.values <= 0 | p.values > 1)) {
    stop("p.values must be in (0, 1].")
  }

  K <- length(p.values)

  p.sorted <- sort(p.values)
  input_names <- names(p.sorted)

  rtp_pvalues <- sapply(1:K, function(j) rtp::p.rtp(j, p.sorted))

  cutpoint <- which.min(rtp_pvalues)
  q_obs <- min(rtp_pvalues)

  log_q <- -log10(q_obs)
  log_q <- pmin(log_q, 15)

  UFI_base <- UFI_database[UFI_database$K == K, ]

  if (nrow(UFI_base) == 0) {
    stop("No UFI table available for K = ", K, ".")
  }

  fit <- stats::lm(ARTP_is ~ stats::poly(q, 3), data = UFI_base)

  log_p <- stats::predict(fit, newdata = data.frame(q = log_q))

  p.val <- min(10^(-log_p), 1)

  list(
    p.value = as.numeric(p.val),
    cutpoint = cutpoint,
    observed_stat = q_obs,
    rtp_pvalues = rtp_pvalues
  )
}




