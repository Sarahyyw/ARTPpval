#' ARTP p-value by cross-entropy importance sampling
#'
#' Estimates an ARTP p-value using cross-entropy importance sampling.
#'
#' @param p.values Numeric vector of observed p-values.
#' @param ro Quantile level for cross-entropy updating. Default is 0.01.
#' @param theta Initial proposal parameter. Default is 1.
#' @param N Number of ISCE samples. Default is 10000.
#' @param idx Decay parameter for mixture probabilities. Default is 1.
#'
#' @return A list with the ARTP p-value, selected cutpoint, observed ARTP statistic,
#' and all RTP p-values across cutpoints.
#'
#' @export
#' @importFrom stats quantile optimize dnorm pnorm rbinom rnorm
#'
#'
artp_isce <- function(p.values, ro = 0.01, theta = 1, N = 10^4, idx = 1) {

  if (!is.numeric(p.values)) {
    stop("p.values must be numeric.")
  }

  if (any(is.na(p.values))) {
    stop("p.values cannot contain NA values.")
  }

  if (any(p.values < 0 | p.values > 1)) {
    stop("p.values must be between 0 and 1.")
  }

  K <- length(p.values)

  rtp_pvals <- sapply(1:K, function(j) rtp::p.rtp(j, p.values))
  cutpoint <- which.min(rtp_pvals)
  q_obs <- min(rtp_pvals)

  ARTP.stat <- function(p.values) {
    K <- length(p.values)
    min(sapply(1:K, function(j) rtp::p.rtp(j, p.values)))
  }

  Samp.mixed.grad <- function(K, theta, N = 1000, idx = 1) {

    grad <- 1 / ((1:K)^idx + 1)

    mix.prop <- matrix(
      rbinom(K * N, size = 1, prob = grad),
      nrow = K
    )
    mix.prop <- t(mix.prop)

    x <- matrix(rnorm(K * N, 0, theta), nrow = N) ^ mix.prop +
      matrix(rnorm(K * N), nrow = N) ^ (1 - mix.prop) - 1

    p.mat <- 2 * pnorm(abs(x), lower.tail = FALSE)
    stat.val <- apply(p.mat, 1, ARTP.stat)

    log_g_norm <- rowSums(log(dnorm(x)))
    log_g_mix <- rowSums(
      log(
        t(t(dnorm(x, 0, theta)) * grad) +
          t(t(dnorm(x)) * (1 - grad))
      )
    )

    w <- exp(log_g_norm - log_g_mix)

    list(stat.val = stat.val, weight = w, x = x)
  }

  mixed.par.update.grad <- function(object, ro, K, idx = 1) {

    stat.val <- object[[1]]
    gamma <- quantile(stat.val, ro)
    weight <- object[[2]]
    x.mat <- t(object[[3]])
    B <- length(stat.val)
    grad <- 1 / ((1:K)^idx + 1)

    opt <- optimize(
      function(t) {
        1 / B * sum(
          weight[stat.val <= gamma] *
            apply(
              as.matrix(x.mat[, stat.val <= gamma]),
              2,
              function(x) {
                sum(log(grad * dnorm(x, 0, t) + (1 - grad) * dnorm(x, 0, 1)))
              }
            )
        )
      },
      interval = c(0, 10),
      maximum = TRUE
    )

    list(gamma = gamma, par = opt$maximum)
  }

  gamma <- Inf
  iter <- 0
  p.val <- NA

  while (gamma > q_obs) {
    iter <- iter + 1

    obj <- Samp.mixed.grad(K = K, theta = theta, N = N, idx = idx)
    stat.val <- obj[[1]]
    gamma <- quantile(stat.val, ro)

    if (iter > 10) {
      p.val <- NA
      break
    }

    if (gamma > q_obs) {
      par <- mixed.par.update.grad(object = obj, ro = ro, K = K, idx = idx)
      theta <- par$par
      gamma <- par$gamma
    }
  }

  message(paste0("t=", iter, ", stop"))

  if (iter <= 10) {
    obj.final <- Samp.mixed.grad(K = K, theta = theta, N = N, idx = idx)
    stat.val.final <- obj.final[[1]]
    weight.final <- obj.final[[2]]

    p.val <- 1 / N * sum(weight.final[stat.val.final <= q_obs])
  }

  list(
    p.value = p.val,
    cutpoint = cutpoint,
    observed_stat = q_obs,
    rtp_pvalues = rtp_pvals
  )
}








