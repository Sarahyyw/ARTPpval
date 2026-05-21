
library(tidyverse)
library(ggplot2)
library(rtp)
library(tidyverse)
library(rlist)
library(parallel)
library(matrixStats)

##### rtp function #####
RTP.stat = function(p.values, i) {
  K=length(p.values)
  v.stat.pval=sapply(i ,function(i) p.rtp(i, p.values))
  return(stat=v.stat.pval)
}

# RTP.stat(c(0.01, 0.02, 0.5, 0.7, 0.9, 0.9), 3)

###### IS for Beta, Exponential and Gaussian #####
ImportSamp=function(k,lambda,family,mu=1,J) {
  if (family=="Beta") {
    x = sapply(lambda, function(x) rbeta(1,shape1 = x,shape2 = 1))
    x_sort = sort(x,decreasing = F)
    stat.val = sum(-log(x_sort[1:J]))
    w = 1/prod(lambda*x^(lambda-1))
    w.ext = w*log(x)
  }
  
  if (family=="Exponential") {
    x = sapply(lambda, function(x) rexp(1,rate =1/x))
    p = pexp(x,rate=1/mu,lower.tail = FALSE)
    p_sort = sort(p, decreasing = F)
    stat.val = sum(-log(p_sort[1:J]))
    w = prod(lambda)/mu^k*exp(sum((1/lambda-1/mu)*x))
    w.ext = w*x
  }
  
  if (family=="Gaussian") {
    x = sapply(lambda, function(x) rnorm(1,mean = 0,sd=sqrt(x)))
    p = 2*pnorm(abs(x),mean=0,sd=sqrt(mu),lower.tail = FALSE)
    p_sort = sort(p, decreasing = F)
    stat.val = sum(-log(p_sort[1:J]))
    w = prod(sqrt(lambda))/mu^(k/2)*exp(sum((1/lambda-1/mu)*x^2/2))
    w.ext = w*x^2
  }
  
  return(list(stat.val=stat.val, w=w, w.ext=w.ext))
}


par.update=function(obj, ro, family) {
  stat.val = sapply(obj,"[[",1)
  gamma = quantile(stat.val,1-ro)
  w = sapply(obj,"[[",2)
  w.ext = sapply(obj,"[[",3)
  
  if (family=="Beta") {
    lambda.num = sum(w[stat.val>=gamma])
    lambda_den = apply(w.ext,1,function(x) sum(x[stat.val>=gamma]))
    lambda = -lambda.num/lambda_den
  }
  
  if (family=="Exponential") {
    lambda.num=apply(w.ext,1,function(x) sum(x[stat.val>=gamma]))
    lambda_den=sum(w[stat.val>=gamma])
    lambda=lambda.num/lambda_den
  }
  
  if (family=="Gaussian") {
    lambda.num=apply(w.ext,1,function(x) sum(x[stat.val>=gamma]))
    lambda_den=sum(w[stat.val>=gamma])
    lambda=lambda.num/lambda_den
  }
  
  return(list(gamma=gamma,lambda=lambda))
}


CE=function(p.val,ro,N=10^5,family,mu=1,J) {
  k=length(p.val)
  p_sort = sort(p.val, decreasing = F)
  q_obs= sum(-log(p_sort[1:J]))
  lambda = rep(mu,k)
  gamma = -Inf
  t = 0
  
  while (gamma<q_obs) {
    t=t+1
    obj = lapply(1:N,function(i) ImportSamp(k=k,lambda=lambda,family=family,mu=mu, J=J))
    stat.val = sapply(obj,"[[",1)
    gamma = quantile(stat.val,1-ro)
    
    if (gamma<q_obs) {
      par = par.update(obj = obj,ro=ro,family=family)
      lambda = par$lambda
      gamma = par$gamma
      #message(paste0(gamma,", ",lambda))
    }
  }
  message("stop")
  message(paste0(t))
  
  obj.final = lapply(1:N,function(i) ImportSamp(k,lambda=lambda,family=family, mu=mu, J=J))
  stat.val.final=sapply(obj.final, "[[",1)
  w.final=sapply(obj.final, "[[",2)
  p.val=1/N*sum(w.final[stat.val.final>=q_obs])
  return(list(p.val=p.val,iter=t,lambda=lambda))
}


########## IS Gaussian Mixture grad #############

Samp.mixed.grad=function(K=k,theta,N=1000,idx=1,J) {
  
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


CE.mixed.grad(p.val=c(0.01, 0.02, 0.5, 0.7, 0.9, 0.9), J=3)
RTP.stat(p.val = c(0.01, 0.02, 0.5, 0.7, 0.9, 0.9), 3)




##### RTP signal_delta_matching target(gamma) #####
#' @param tol     Tolerance for uniroot
#' @param maxit   Max iterations for uniroot
#' @param s_lower Lower bound for s (close to 0 -> very strong signals)
#' @param s_upper Upper bound for s (close to 1 -> very weak signals)
#' @return list with s_signal, achieved_RTP_p, pvalues_example, etc.
#' 
rtp_match_signal <- function(K = 50,
                             S = 5,
                             target = 1e-3,
                             noise_p = NULL,
                             seed = 990123,
                             tol = 1e-10,
                             maxit = 100,
                             delta_lower = 1e-16,
                             delta_upper = 0.999999) {
  # set.seed(seed)
  noise_p <- runif(K - S)
  
  # objective: RTP(m, [s,..,s, noise]) - target
  f <- function(delta) {
    # guard for boundary issues
    if (!is.finite(delta) || delta <= 0 || delta >= 1) return(NA_real_)
    pvals <- c(rep(delta, S), noise_p)
    p.rtp(S, pvals) - target
  }
  # bracket check
  f_lo <- f(delta_lower); f_hi <- f(delta_upper)
  
  if (is.na(f_lo) || is.na(f_hi)) {stop("p.rtp returned NA at the search bounds; adjust s_lower/s_upper.")}
  if (f_lo > 0 && f_hi > 0) {stop("Target too small to reach (even at very strong signals). Try smaller s_lower or revise target.")}
  if (f_lo < 0 && f_hi < 0) {stop("Target too large to reach (even at very weak signals). Try larger s_upper or revise target.")}
  # solve
  root <- uniroot(f, lower = delta_lower, upper = delta_upper, tol = tol, maxiter = maxit)
  delta_hat <- root$root
  # achieved <- p.rtp(m, c(rep(delta_hat, m), noise_p))
  # pvals_example <- c(rep(delta_hat, m), noise_p)
  return(delta_hat)
}





##### generate input p-value #####

beta =  seq(from = 0, to = 1, by = 0.1)
gamma = c(10^(-3), 10^(-5), 10^(-7)) ## target p-value
## define input p-values - for several number of num_sig
input_generate <- function(n, beta, gamma) {
  s = ceiling(n^(beta))
  p_input_beta = c()
  signal_sum = data.frame(num_sig = rep(s, each = 3) , target = rep(gamma, 11), delta = rep(NA,33))
  
  for (i in 1:length(s)) {
    s_i = s[i]
    p_noise = runif(n-s_i, 0, 1)
    
    p_input_gamma = c()
    for (j in 1:length(gamma)) {
      p_target = gamma[j]
      
      # t = qchisq(p = 1 - p_target, df = 2*n)
      # delta = exp(-(t + 2*sum(log(p_noise)))/(2*s_i))
      delta = rtp_match_signal(K=n, S=s_i, target=p_target)
      
      p_input_j = c(rep(delta, s_i), p_noise)
      signal_sum$delta[j+((i-1)*3)] = delta
      p_input_gamma[[j]] = p_input_j
      
    }
    p_input_beta[[i]] = p_input_gamma
  }
  final = list(signal_sum = signal_sum, p_input = p_input_beta)
  return(final)
}







##### run for rtp #####
## process one row in the res data frame, the idx is just a parameter in the gaussian mixture, mostly just use idx=1

process_row <- function(i, res, N, p_input=p_input, n_input=n_input, sig_unique=sig_unique, gamma=gamma) {
  num_sig = res$num_sig[i]
  beta_i = log(num_sig, base=n_input)
  input_i = p_input$p_input[[which(sig_unique ==  num_sig)]]
  gamma_cur = res$gamma[i]
  j = which(gamma == gamma_cur)
  input_i_j = input_i[[j]]
  
  if (beta_i < 0.5) {idx = 1
  } else if (beta_i >= 0.5) {idx=0.3}
  
  res$idx[i] = idx; res$N[i] = N
  res$p_rtp[i] = RTP.stat(input_i_j, res$check_point[i])   
  rtp_cur <- mclapply(1:10, function(l) {
    CE.mixed.grad(p.val = input_i_j, ro = 0.01, J = res$check_point[i], theta = 1, N = N, idx = idx)}, mc.cores = 5)
  res[i, grep("^p_is", colnames(res))] <- unlist(rtp_cur)
  message("num_sig=", num_sig)
  return(res[i, ])
}




run_for_RTP <- function(p_input = p_input_50, N=1e2, n_input=50) {
  
  n_max = n_input
  N = N
  sig_sum = p_input$signal_sum
  sig = sig_sum$num_sig
  sig_unique = unique(sig)
  check_point = pmin(c(rbind(sig, sig * 2)), n_max)
  family = family
  target = sig_sum$target
  delta = sig_sum$delta
  
  
  ### for each num_sig and each delta, check j and 2*j for rtp
  res = data.frame(num_sig = rep(sig_sum$num_sig, each=2), gamma = rep(target, each=2), 
                   check_point = check_point, idx = NA, N = NA,
                   p_rtp = NA, p_is_1 = NA, p_is_2 = NA, p_is_3 = NA, p_is_4 = NA, p_is_5 = NA, p_is_6 = NA, 
                   p_is_7 = NA, p_is_8 = NA, p_is_9 = NA, p_is_10 = NA) 
  
  res_sub <- mclapply(1:nrow(res), function(i)  {
    process_row(i=i,res=res, N=N, p_input=p_input, n_input=n_input, sig_unique=sig_unique, gamma=gamma)
    }, mc.cores=4)
  
  res <- do.call(rbind, res_sub)
  return(res)
}



######## Benchmark & tau bound #####
benchmark_func = function(res) {
  
  cols <- paste0("p_is_", 1:10)
  p_is_mean <- rowMeans(res[, cols], na.rm = TRUE)
  res$p_is_mean = p_is_mean
  p_is_var <- rowVars(as.matrix(res[, cols]), na.rm = TRUE)

  first_term_bench = abs(res$p_rtp - p_is_mean)/res$p_rtp
  second_term_bench = sqrt(100*p_is_var/(res$p_rtp^2))
  res$benchmark = first_term_bench + second_term_bench
  res = res %>% dplyr::select(p_rtp, p_is_mean, benchmark, N)
  return(res)
}

benchmark_tau <- function(p_m) {
  if (p_m <= 1e-7) {
    return(3.0)
  } else if (p_m <= 1e-5 & p_m > 1e-7) {
    return(1.0)
  } else if (p_m < 1e-3 & p_m > 1e-5) {
    return(0.3)
  } else {
    return(0.10)
  }
}


