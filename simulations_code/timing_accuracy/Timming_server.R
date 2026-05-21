library(tidyverse)
library(ggplot2)
library(rtp)
library(tidyverse)
library(rlist)
library(parallel)

##### IS Gaussian Mixture grad #####

Samp.mixed.grad = function(K,theta,N=1000,idx=1) {
  
  grad=1/((1:K)^idx+1)
  mix.prop=matrix(rbinom(K*N,size=1,prob=grad),nrow=K)
  mix.prop=t(mix.prop)
  x=matrix(rnorm(K*N,0,theta),nrow = N) ^ {(mix.prop)} + matrix(rnorm(K*N),nrow = N) ^ {1-(mix.prop)} - 1
  
  p.mat = 2*pnorm(abs(x),lower.tail = FALSE)
  stat.val = apply(p.mat, 1, RTP.vse)
  
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
  #message(paste0("t=", t, ", stop"))
  
  if (t <= 10) {
    objs = lapply(1, function(b) Samp.mixed.grad(K, theta = theta, N = N, idx = idx))
    stat.val.final = as.vector(sapply(objs, "[[", 1))
    weight.final = as.vector(sapply(objs, "[[", 2))
    p.val = 1/N * sum(weight.final[stat.val.final <= q_obs])
  }
  
  return(p.val)
}

# ##### first layer MC, Second layer MC ########
# RTP.MC=function(p.values, N=1e3) {
#   K=length(p.values)
#   p.values=sort(p.values,decreasing = F)
#   q.values=cumsum(-log(p.values))
#   p.mat=matrix(runif(K*N),nrow = N)
#   q.mat=t(apply(p.mat,1,function(x) {
#     cumsum(-log(sort(x,decreasing = F)))
#   }))
#   #print(dim(q.mat))
#   q.mat=t(apply(q.mat,1,function(x) x >= q.values))
#   stat=min((apply(q.mat,2, function(x) (sum(x)+1)/(N+1))))
#   return(stat=stat)
# }
# 
# ARTP.MC.MC = function(K, N=10^3, q) {
#   stat.vals=lapply(1:N, function(n) {
#     p0 = runif(K)
#     stat.val = RTP.MC(p0,N=N)
#     return(stat.val)
#   })
#   stat.vals=unlist(stat.vals)
#   p=mean(stat.vals <= q)
#   return(p)
# } 

##### r.rtp for inner layer, MC for outer layer #######
RTP.vse = function(p.values) {
  K = length(p.values)
  v.stat.pval = sapply(1:K,function(i) p.rtp(i, p.values))
  v.AFp = min(v.stat.pval)
  return(stat=v.AFp)
}

ARTP.vse.MC = function(K, N=10^4, q) {
  stat.vals=lapply(1:N, function(n) {
    p0 = runif(K)
    stat.val = RTP.vse(p0)
    return(stat.val)
  })
  stat.vals=unlist(stat.vals)
  p=mean(stat.vals <= q)
  return(p)
} 


####### r.rtp for inner layer, IS for outer layer #######
ARTP.vse.IS <- CE.mixed.grad


####### r.rtp for inner layer, UFI for outer layer #######

ARTP.vse.UFI = function(q, k) {
  load("/Users/sarah/Desktop/ARTP/result/UFI_database.RData")
  UFI_database = res
  UFI = UFI_database %>% filter(K == K)
  polyd3 <- lm(ARTP_is ~ poly(q, 3), data = UFI)
  q_new = -log10(q)
  p_UFI <- 10^(-predict(polyd3, newdata = data.frame(q=q_new)))
  p_UFI = ifelse(p_UFI>1, 1, p_UFI)
  return(p_UFI)
}


#### time running #####
time_individual_10 = function(r, res) {
  res_r = res[r,]
  q=res_r$q; K=res_r$K; N=res_r$N; method=res_r$method; idx=res_r$idx
  
  # if (method == "MC_MC") {
  #   start_time=Sys.time()
  #   p_10 = mclapply(1:10, function(i) {ARTP.MC.MC(K, N=N, q)},mc.cores=10)
  #   end_time <- Sys.time()
  #   t=as.numeric(difftime(end_time, start_time, units = "secs"))
  #   as.numeric(difftime(end_time, start_time, units = "secs"))
  # }
  
  if (method == "vse_MC") {
    start_time=Sys.time()
    p_10 = mclapply(1:10, function(i) {ARTP.vse.MC(K, N=N, q)},mc.cores=10)
    end_time <- Sys.time()
    t=as.numeric(difftime(end_time, start_time, units = "secs"))
  }
  
  if (method == "vse_IS") {
    start_time=Sys.time()
    p_10 = mclapply(1:10, function(i) {ARTP.vse.IS(q_obs=q,K=K,ro=0.01,theta=1,N=N,idx=idx)}, mc.cores=10)
    end_time <- Sys.time()
    t=as.numeric(difftime(end_time, start_time, units = "secs"))
  }
  
  if (method == "vse_UFI") {
    start_time=Sys.time()
    p_10 = mclapply(1:10, function(i) {ARTP.vse.UFI(q=q,K=K)},mc.cores=10)
    end_time <- Sys.time()
    t=as.numeric(difftime(end_time, start_time, units = "secs"))
  }
  
  message(paste0("q=", q, ", K=", K, ", method=", method))
  
  p_10 = unlist(p_10)
  res_with_time = (list(time=t, p_10=p_10))
  return(res_with_time)
}


time_individual_1 = function(r, res) {
  res_r = res[r,]
  q=res_r$q; K=res_r$K; N=res_r$N; method=res_r$method; idx=res_r$idx
  
  # if (method == "MC_MC") {
  #   start_time=Sys.time()
  #   p_1 = ARTP.MC.MC(K, N=N, q)
  #   end_time <- Sys.time()
  #   t=as.numeric(difftime(end_time, start_time, units = "secs"))
  #   as.numeric(difftime(end_time, start_time, units = "secs"))
  # }
  
  if (method == "vse_MC") {
    start_time=Sys.time()
    p_1 = ARTP.vse.MC(K, N=N, q)
    end_time <- Sys.time()
    t=as.numeric(difftime(end_time, start_time, units = "secs"))
  }
  
  if (method == "vse_IS") {
    start_time=Sys.time()
    p_1 = ARTP.vse.IS(q_obs=q,K=K,ro=0.01,theta=1,N=N,idx=idx)
    end_time <- Sys.time()
    t=as.numeric(difftime(end_time, start_time, units = "secs"))
  }
  
  if (method == "vse_UFI") {
    start_time=Sys.time()
    p_1 = ARTP.vse.UFI(q=q, K=K)
    end_time <- Sys.time()
    t=as.numeric(difftime(end_time, start_time, units = "secs"))
  }
  message(paste0("q=", q, ", K=", K, ", method=", method))
  
  p_1 = p_1
  res_with_time = (list(time=t, p_1=p_1))
  return(res_with_time)
}


time_ARTP_10 = function(q_input=q_input,K=K,method=method,idx=idx,N=N) {
  
  res = data.frame(q=q_input, K=K, method=method, idx=idx, N=N, time=NA,
                   p_1=NA,p_2=NA,p_3=NA,p_4=NA,p_5=NA,p_6=NA,p_7=NA,p_8=NA,p_9=NA,p_10=NA)
  
  res_list = mclapply(1:nrow(res), function(r) {time_individual_10(r, res)}, mc.cores=1)
  res$time = unlist(lapply(res_list, function(x) x$time))
  p_10_all = lapply(res_list, function(x) x$p_10)
  for (i in 1:nrow(res)) {res[i, paste0("p_", 1:10)] <- p_10_all[[i]]}
  return(res)
}


time_ARTP_1 = function(q_input=q_input,K=K,method=method,idx=idx,N=N) {
  
  res = data.frame(q=q_input, K=K, method=method, idx=idx, N=N, time=NA,p_1=NA)
  res_list = mclapply(1:nrow(res), function(r) {time_individual_1(r, res)}, mc.cores=15)
  res$time = unlist(lapply(res_list, function(x) x$time))
  res$p_1 = unlist(lapply(res_list, function(x) x$p_1))
  return(res)
}


