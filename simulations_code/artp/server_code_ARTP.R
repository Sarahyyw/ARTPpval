library(tidyverse)
library(ggplot2)
library(rtp)
library(tidyverse)
library(rlist)
library(parallel)
# core.num = 30

##### rtp function #####
RTP.stat = function(p.values, i) {
  K=length(p.values)
  v.stat.pval = sapply(i ,function(i) p.rtp(i, p.values))
  return(stat=v.stat.pval)
}
## the p.rtp function itself will sort the data first, no need to sort it in IS functions

ARTP.stat = function(p.values) {
  K = length(p.values)
  v.stat.pval = sapply(1:K,function(i) p.rtp(i, p.values))
  v.AFp = min(v.stat.pval)
  return(stat=v.AFp)
}

##### IS for Beta, Exponential and Gaussian #####
ImportSamp=function(K,lambda,family,mu=1) {
  if (family=="Beta") {
    x = sapply(lambda, function(x) rbeta(1,shape1 = x,shape2 = 1))
    w = 1/prod(lambda*x^(lambda-1)) 
    w.ext = w*log(x) 
    stat.val = ARTP.stat(p.values=x)
  }
  
  if (family=="Exponential") {
    x = sapply(lambda, function(x) rexp(1, rate =1/x))
    p = pexp(x,rate=1/mu,lower.tail = FALSE)
    w = prod(lambda)/mu^K*exp(sum((1/lambda-1/mu)*x))
    w.ext = w*x
    stat.val =  ARTP.stat(p.values=p)
  }
  
  if (family=="Gaussian") {
    x = sapply(lambda, function(x) rnorm(1,mean = 0,sd=sqrt(x)))
    p = 2*pnorm(abs(x),mean=0,sd=sqrt(mu),lower.tail = FALSE)
    w = prod(sqrt(lambda))/mu^(K/2)*exp(sum((1/lambda-1/mu)*x^2/2))
    w.ext = w*x^2
    stat.val =  ARTP.stat(p.values=p)
  }
  
  return(list(stat.val=stat.val, w=w, w.ext=w.ext))
}


par.update=function(obj, ro, family) {
  stat.val = sapply(obj,"[[",1)
  gamma = quantile(stat.val, ro)
  w = sapply(obj,"[[",2)
  w.ext = sapply(obj,"[[",3)
  
  if (family=="Beta") {
    lambda.num = sum(w[stat.val <= gamma])
    lambda_den = apply(w.ext,1,function(x) sum(x[stat.val <= gamma]))
    lambda = -lambda.num/lambda_den
  }
  
  if (family=="Exponential") {
    lambda.num=apply(w.ext,1,function(x) sum(x[stat.val <= gamma]))
    lambda_den=sum(w[stat.val <= gamma])
    lambda=lambda.num/lambda_den
  }
  
  if (family=="Gaussian") {
    lambda.num=apply(w.ext,1,function(x) sum(x[stat.val <= gamma]))
    lambda_den=sum(w[stat.val <= gamma])
    lambda=lambda.num/lambda_den
  }
  
  return(list(gamma=gamma,lambda=lambda))
}


CE=function(q_obs,K,ro,N=10^5,family,mu=1) {
  
  lambda = rep(mu,K)
  gamma = Inf
  t = 0
  
  while (gamma > q_obs) {
    t=t+1
    obj = mclapply(1:N,function(i) ImportSamp(lambda=lambda,K,family=family,mu=mu), mc.cores=core.num)
    stat.val = sapply(obj,"[[",1)
    gamma = quantile(stat.val, ro)
    
    if (gamma > q_obs) {
      par = par.update(obj = obj,ro=ro,family=family)
      lambda = par$lambda
      gamma = par$gamma
      # message(paste0(gamma,", ",lambda))
    }
  }
  message("stop")
  
  obj.final = mclapply(1:N,function(i) ImportSamp(lambda=lambda,family=family,K,mu=mu), mc.cores=core.num)
  stat.val.final=sapply(obj.final, "[[",1)
  w.final=sapply(obj.final, "[[",2)
  p.val = 1/N*sum(w.final[stat.val.final <= q_obs])
  
  return(list(p.val=p.val,iter=t,lambda=lambda))
}

##### IS Gaussian Mixture grad #####

Samp.mixed.grad = function(K,theta,N=1000,idx=1) {
  
  grad=1/((1:K)^idx+1)
  mix.prop=matrix(rbinom(K*N,size=1,prob=grad),nrow=K)
  mix.prop=t(mix.prop)
  x=matrix(rnorm(K*N,0,theta),nrow = N) ^ {(mix.prop)} + matrix(rnorm(K*N),nrow = N) ^ {1-(mix.prop)} - 1
  
  p.mat = 2*pnorm(abs(x),lower.tail = FALSE)
  stat.val = apply(p.mat, 1, ARTP.stat)
  
  log_g_norm <- rowSums(log(dnorm(x))) # numerator part
  log_g_mix <- rowSums(log(t(t(dnorm(x,0,theta))*grad)+t(t(dnorm(x))*(1-grad))))
  w = exp(log_g_norm-log_g_mix)
  #s=rowsums(x[,grad==1]^2)
  return(list(stat.val=stat.val, weight=w, x=x))
}


mixed.par.update.grad=function(object,ro,K,idx=1) {
  
  stat.val=object[[1]]
  gamma=quantile(stat.val, ro)
  weight=object[[2]]
  x.mat=t(object[[3]])
  B=length(stat.val)
  grad=1/((1:K)^idx+1)
  
  opt=optimize(function(t) 1/B*sum(weight[stat.val <= gamma]*
                                     apply(as.matrix(x.mat[,stat.val <= gamma]), 2, function(x) sum(log(grad*dnorm(x,0,t)+(1-grad)*dnorm(x,0,1))))), 
               interval = c(0,10), maximum = TRUE)
  
  par=opt$maximum
  return(list(gamma=gamma,par=par))
}


CE.mixed.grad = function(q_obs,K, ro=0.01, theta=1, N=10^4, idx=1) {
  
  theta = theta; idx = idx; N = N
  gamma = Inf
  t = 0
  
  while (gamma > q_obs) {
    t=t+1
    obj=Samp.mixed.grad(K, theta=theta, N=N, idx=idx)
    stat.val=obj[[1]]
    gamma = quantile(stat.val, ro)
    
    if (t > 10) {
      p.val = NA  
      break
    }
    
    if (gamma > q_obs) {
      par = mixed.par.update.grad(object = obj, ro=ro, K, idx = idx)
      theta=par$par
      gamma=par$gamma
      #message(paste0(gamma,", ",theta))
    }
  }
  message(paste0("t=", t, ", stop"))
  
  if (t <= 10) {
    objs = lapply(1, function(b) Samp.mixed.grad(K, theta = theta, N = N, idx = idx))
    stat.val.final = as.vector(sapply(objs, "[[", 1))
    weight.final = as.vector(sapply(objs, "[[", 2))
    p.val = 1/N * sum(weight.final[stat.val.final <= q_obs])
  }
  
  return(p.val)
}


##### IS Gaussian Mixture prop #####

Samp.mixed.prop=function(K,theta,N1=1000,k0,k1,thre,prop=0.2) {
  
  if (prop==-Inf) {grad=rep(1,K)} else {
    grad=c(rep(0,ceiling(K^(prop))),rep(1,floor(K-K^(prop))))}
  
  #grad=1/((1:K)^(0.1)+1)
  mix.prop=matrix(rbinom(K*N,size=1,prob=grad),nrow=K)
  mix.prop=t(mix.prop)
  x=matrix(rnorm(K*N,0,theta),nrow = N)^{(mix.prop)}+matrix(rnorm(K*N),nrow = N)^{1-(mix.prop)}-1
  
  p.mat = 2*pnorm(abs(x),lower.tail = FALSE)
  stat.val = apply(p.mat, 1, ARTP.stat)
  
  log_g_norm <- rowSums(log(dnorm(x))) # numerator part
  log_g_mix <- rowSums(log(t(t(dnorm(x,0,theta))*grad)+t(t(dnorm(x))*(1-grad))))
  w=exp(log_g_norm-log_g_mix)
  
  return(list(stat.val=stat.val,weight=w,x=x))
}


mixed.par.update.prop=function(object,ro,K,prop) {
  
  stat.val=object[[1]]
  gamma=quantile(stat.val, ro)
  weight=object[[2]]
  x.mat=t(object[[3]])
  B=length(stat.val)
  
  if (prop==-Inf) {grad=rep(1,K)} else {
    grad=c(rep(0,ceiling(K^(prop))),rep(1,floor(K-K^(prop))))}
  
  opt=optimize(function(t) 1/B*sum(weight[stat.val <= gamma]*
                                     apply(as.matrix(x.mat[,stat.val <= gamma]), 2, function(x) sum(log(grad*dnorm(x,0,t)+(1-grad)*dnorm(x,0,1))))), 
               interval = c(0,5), maximum = TRUE)
  
  par=opt$maximum
  return(list(gamma=gamma,par=par))
}


CE.mixed.prop=function(q_obs, K, ro, N, prop=0.2, theta=1) {
  
  gamma=Inf
  t=0
  
  while (gamma > q_obs) {
    t = t+1
    obj = Samp.mixed.prop(K=K, theta=theta, N=N, prop = prop)
    stat.val=obj[[1]]
    gamma = quantile(stat.val, ro)
    
    if (t > 10) {
      p.val = NA  
      break
    }
    
    if (gamma > q_obs) {
      par = mixed.par.update.prop(object = obj, ro=ro, K, prop=prop)
      theta = par$par
      gamma=par$gamma
      # message(paste0(gamma,", ",theta))
    }
  }
  print("stop")
  
  if (t <= 10) {
    objs = lapply(1, function(b) Samp.mixed.prop(K, theta = theta, N = N, prop = prop))
    stat.val.final = as.vector(sapply(objs, "[[", 1))
    weight.final = as.vector(sapply(objs, "[[", 2))
    p.val = 1/N * sum(weight.final[stat.val.final <= q_obs])
  }
  
  return(p.val)
}


##### running function ######
process_row <- function(r, res, K, N, idx, prop, family) {
  q <- res$q[r]
  if (family == "Gaussian mixture grad") {
    artp_cur <- CE.mixed.grad(q_obs = q, K = K, ro = 0.01, theta = 1, N = N, idx = idx)
  }
  return(artp_cur)
}

# process_row(r,res=res, K=K,N=N,idx=idx,prop=prop, family=family)
##### running function ######
run_for_p_value_ARTP <- function(q_input=q_input, K=50, N=1e3, family= "Exponential", idx=idx, prop=prop) {
  
  q_input = q_input
  family = family
  
  res = data.frame(q=q_input, K=K, family=family, idx=idx, N=N, prop = prop, ARTP_is = NA)
  
  if (family == "Gaussian mixture grad" | family == "Gaussian mixture prop") {
    res_sub <- mclapply(1:nrow(res), function(r)  process_row(r,res=res, K=K,N=N,idx=idx,prop=prop, family=family),  mc.cores=15)
    res$ARTP_is <- unlist(res_sub)
  }
  
  return(res)
}

   




