
#' @importFrom stats quantile optimize dnorm pnorm rbinom rnorm
NULL

########## IS Gaussian Mixture grad #############

Samp.mixed.grad=function(K,theta,N=1000,idx=1,J) {

  grad=1/((1:K)^idx+1)
  mix.prop=matrix(rbinom(K*N,size=1,prob=grad),nrow=K)
  mix.prop=t(mix.prop)
  x=matrix(rnorm(K*N,0,theta),nrow = N) ^ {(mix.prop)} + matrix(rnorm(K*N),nrow = N) ^ {1-(mix.prop)} - 1

  p.mat = 2*pnorm(abs(x),lower.tail = FALSE)
  p_sort.mat = Rfast::rowSort(p.mat, descending = F)
  stat.val = rowSums(as.matrix(-log(p_sort.mat[,1:J])))

  log_g_norm <- rowSums(log(dnorm(x))) # numerator part
  log_g_mix <- rowSums(log(t(t(dnorm(x,0,theta))*grad)+t(t(dnorm(x))*(1-grad))))
  w=exp(log_g_norm-log_g_mix)
  #s=rowsums(x[,grad==1]^2)
  return(list(stat.val=stat.val, weight=w, x=x))
}


mixed.par.update.grad=function(object,ro,K,idx=1) {

  stat.val=object[[1]]
  gamma=quantile(stat.val,1-ro)
  weight=object[[2]]
  x.mat=t(object[[3]])
  B=length(stat.val)
  grad=1/((1:K)^idx+1)

  opt=optimize(function(t) 1/B*sum(weight[stat.val>= gamma]*
                                     apply(as.matrix(x.mat[,stat.val>=gamma]), 2, function(x) sum(log(grad*dnorm(x,0,t)+(1-grad)*dnorm(x,0,1))))),
               interval = c(0,10), maximum = TRUE)

  par=opt$maximum
  return(list(gamma=gamma,par=par))
}


#' RTP p-value by cross-entropy importance sampling
#'
#' Estimates an RTP p-value using cross-entropy importance sampling.
#'
#' @param p.val Numeric vector of p-values.
#' @param ro Quantile level for cross-entropy updating.
#' @param J Number of smallest p-values used in the RTP statistic.
#' @param theta Initial proposal parameter.
#' @param N Number of ISCE sampling.
#' @param idx Decay parameter for mixture probabilities.
#'
#' @return Estimated RTP p-value.
#'
#' @export

CE.mixed.grad = function(p.val,ro=0.01, J, theta=1, N=10^4, idx=1) {

  k=length(p.val)
  p_sort = sort(p.val, decreasing = F)
  q_obs= sum(-log(p_sort[1:J]))
  theta = theta; idx = idx; N = N
  gamma = -Inf
  t = 0

  while (gamma<q_obs) {
    t=t+1
    obj=Samp.mixed.grad(K=k, theta=theta, N=N, idx=idx, J = J)
    stat.val=obj[[1]]
    gamma = quantile(stat.val,1-ro)

    if (t > 20) {
      p.val = NA
      break
    }

    if (gamma<q_obs) {
      par = mixed.par.update.grad(object = obj, ro=ro, K=k, idx = idx)
      theta=par$par
      gamma=par$gamma
      #message(paste0(gamma,", ",theta))
    }
  }
  message(paste0("t=", t, ", stop"))

  if (t <= 10) {
    objs = lapply(1, function(b) Samp.mixed.grad(K = k, theta = theta, N = N, idx = idx, J = J))
    stat.val.final = as.vector(sapply(objs, "[[", 1))
    weight.final = as.vector(sapply(objs, "[[", 2))
    p.val = 1 / N * sum(weight.final[stat.val.final >= q_obs])
  }

  return(p.val)
}

#' RTP p-value by cross-entropy importance sampling
#'
#' This is an alias for [CE.mixed.grad()].
#'
#' @param p.val Numeric vector of p-values.
#' @param ro Quantile level for cross-entropy updating.
#' @param J Number of smallest p-values used in the RTP statistic.
#' @param theta Initial proposal parameter.
#' @param N Number of ISCE samples.
#' @param idx Decay parameter for mixture probabilities.
#'
#' @return Estimated RTP p-value.
#'
#' @export
rtp_isce <- function(p.val, ro = 0.01, J, theta = 1, N = 10^4, idx = 1) {
  CE.mixed.grad(p.val = p.val, ro = ro, J = J, theta = theta, N = N, idx = idx)
}


# CE.mixed.grad(p.val=c(0.01, 0.02, 0.5, 0.7, 0.9, 0.9), J=3)
# RTP.stat(p.val = c(0.01, 0.02, 0.5, 0.7, 0.9, 0.9), 3)

