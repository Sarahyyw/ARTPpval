## update UFI for COVID up to 1082
load("/Users/sarah/Desktop/ARTP/result/UFI_database.RData")
load("/Users/sarah/Desktop/ARTP/result/ARTP/UFI/ARTP_1080_1088_res_grad_1e4.RData")

ARTP_1080_1088 <- lapply(ARTP_1080_1088_res_grad_1e4, function(x) {return(x$res)})
ARTP_1080_1088 <- do.call(rbind, ARTP_1080_1088)
ARTP_1080_1088 [,c(1,7)] = -log10(ARTP_1080_1088 [,c(1,7)])
UFI_database = rbind(res, ARTP_1080_1088)
# max(UFI_base$K)


############### Geog Region ########################
pval_by_location = readRDS("/Users/sarah/Desktop/ARTP/Real Application/P_individual_Geog_region.rds")
library(rtp)
library(tidyverse)

### getting ARTP 
ARTP.q = function(p.values) {
  K = length(p.values)
  v.stat.pval = sapply(1:K,function(i) p.rtp(i, p.values))
  v.AFp = min(v.stat.pval)
  return(stat=v.AFp)
}

region_ARTP_P = c()
for (i in 1:length(pval_by_location)) {
  region_p=pval_by_location[i][[1]]
  ARTP_p=c()
  for (j in 1:5) {
    region_p_T = region_p[[j]]
    region_q_T = ARTP.q(region_p_T)
    region_q_T = -log10(region_q_T)
    region_q_T = ifelse(region_q_T>15, 15, region_q_T)
    ### get ARTP p-value
    UFI_base = UFI_database %>% filter(K==length(region_p_T))
    polyd3 <- lm(ARTP_is ~ poly(q, 3), data = UFI_base)
    p_T_UFI <- predict(polyd3, newdata = data.frame(q=region_q_T))
    ARTP_p = c(ARTP_p, p_T_UFI)
  }
  region_ARTP_P = c(region_ARTP_P, list(ARTP_p))
}


## Get Fisher
library(poolr)
region_Fisher_P = c()
for (i in 1:length(pval_by_location)) {
  region_p=pval_by_location[i][[1]]
  Fisher_p=c()
  for (j in 1:5) {
    region_p_T = region_p[[j]]
    Fisher_p_T = (fisher(region_p_T))$p
    Fisher_p = c(Fisher_p, Fisher_p_T)
  }
  region_Fisher_P = c(region_Fisher_P, list(Fisher_p))
}


## Get HC
library(HCp)
region_HC_P = c()
for (i in 1:length(pval_by_location)) {
  region_p=pval_by_location[i][[1]]
  HC_p=c()
  for (j in 1:5) {
    region_p_T = region_p[[j]]
    HC_q_T = HCstat(region_p_T, k0=1, thre=FALSE); HC_q_T
    HC_p_T = ufi.p(flibs=HC_flibs, K=length(pval_by_location), q=HC_q_T); HC_p_T
    HC_p_T <- ifelse(is.nan(HC_p_T), 1, HC_p_T); HC_p_T
    HC_p = c(HC_p, HC_p_T)
  }
  region_HC_P = c(region_HC_P, list(HC_p))
}


names(region_ARTP_P) = names(pval_by_location)
region_plot = as.data.frame(do.call(rbind, region_ARTP_P))
colnames(region_plot) = c("T1","T2","T3","T4","T5"); region_plot$region = rownames(region_plot)

region_plot = pivot_longer(
  region_plot,
  cols = starts_with("T"),        
  names_to = "Period",
  values_to = "-log10P"
)

region_plot$`-log10P` = ifelse(region_plot$`-log10P`<0, 0, region_plot$`-log10P`)
region_plot$ARTP = region_plot$`-log10P`
region_plot$Fisher = ifelse(-log10(unlist(region_Fisher_P))>15, 15, -log10(unlist(region_Fisher_P)))
region_plot$HC = ifelse(-log10(unlist(region_HC_P))>15, 15, -log10(unlist(region_HC_P)))
region_plot = region_plot %>% select(-`-log10P`)

region_plot = pivot_longer(region_plot, cols=c(ARTP, Fisher, HC), 
                           names_to = "method", 
                           values_to = "p_value")

Geog_region = region_plot %>% 
  ggplot(aes(x = Period, y = p_value, color = method, shape = method, group = method)) + 
  geom_point(size = 1) + 
  geom_line(alpha = 0.5) + 
  labs(x = "Time Window", y = "-log_10 P-value", title = "Geographic Region") + 
  facet_wrap(~region, nrow = 1) + 
  theme_bw() +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 12),  
    axis.title.x = element_text(size = 10),               
    axis.title.y = element_text(size = 10),               
    legend.position = "none"                           
  )



############### BEA Region ########################
pval_by_location = readRDS("/Users/sarah/Desktop/ARTP/Real Application/P_individual_BEA_region.rds")
library(rtp)
library(tidyverse)


### getting ARTP 
ARTP.q = function(p.values) {
  K = length(p.values)
  v.stat.pval = sapply(1:K,function(i) p.rtp(i, p.values))
  v.AFp = min(v.stat.pval)
  return(stat=v.AFp)
}
region_ARTP_P = c()

for (i in 1:length(pval_by_location)) {
  region_p=pval_by_location[i][[1]]
  ARTP_p=c()
  for (j in 1:5) {
    region_p_T = region_p[[j]]
    region_q_T = ARTP.q(region_p_T)
    region_q_T = -log10(region_q_T)
    region_q_T = ifelse(region_q_T>15, 15, region_q_T)
    k = length(region_p_T)
    
    ### get ARTP p-value
    # ------- UFI ------ #
    UFI_base = UFI_database %>% filter(K==k)
    polyd3 <- lm(ARTP_is ~ poly(q, 3), data = UFI_base)
    p_T_UFI <- predict(polyd3, newdata = data.frame(q=region_q_T))
    ARTP_p = c(ARTP_p, p_T_UFI)
  }
  region_ARTP_P = c(region_ARTP_P, list(ARTP_p))
}


## Get Fisher
library(poolr)
region_Fisher_P = c()
for (i in 1:length(pval_by_location)) {
  region_p=pval_by_location[i][[1]]
  Fisher_p=c()
  for (j in 1:5) {
    region_p_T = region_p[[j]]
    Fisher_p_T = (fisher(region_p_T))$p
    Fisher_p = c(Fisher_p, Fisher_p_T)
  }
  region_Fisher_P = c(region_Fisher_P, list(Fisher_p))
}


## Get HC
library(HCp)
region_HC_P = c()
for (i in 1:length(pval_by_location)) {
  region_p=pval_by_location[i][[1]]
  HC_p=c()
  for (j in 1:5) {
    region_p_T = region_p[[j]]
    HC_q_T = HCstat(region_p_T, k0=1, thre=FALSE); HC_q_T
    HC_p_T = ufi.p(flibs=HC_flibs, K=length(pval_by_location), q=HC_q_T); HC_p_T
    HC_p_T <- ifelse(is.nan(HC_p_T), 1, HC_p_T); HC_p_T
    HC_p = c(HC_p, HC_p_T)
  }
  region_HC_P = c(region_HC_P, list(HC_p))
}


names(region_ARTP_P) = names(pval_by_location)
region_plot = as.data.frame(do.call(rbind, region_ARTP_P))
colnames(region_plot) = c("T1","T2","T3","T4","T5"); region_plot$region = rownames(region_plot)

region_plot = pivot_longer(
  region_plot,
  cols = starts_with("T"),        
  names_to = "Period",
  values_to = "-log10P"
)

region_plot$`-log10P` = ifelse(region_plot$`-log10P`<0, 0, region_plot$`-log10P`)
region_plot$ARTP = region_plot$`-log10P`
region_plot$Fisher = ifelse(-log10(unlist(region_Fisher_P))>15, 15, -log10(unlist(region_Fisher_P)))
region_plot$HC = ifelse(-log10(unlist(region_HC_P))>15, 15, -log10(unlist(region_HC_P)))
region_plot = region_plot %>% select(-`-log10P`)

region_plot = pivot_longer(region_plot, cols=c(ARTP, Fisher, HC), 
                           names_to = "method", 
                           values_to = "p_value")

BEA_region = region_plot %>% 
  ggplot(aes(x = Period, y = p_value, color = method, shape = method, group = method)) + 
  geom_point(size = 1) + 
  geom_line(alpha = 0.5) + 
  labs(x = "Time Window", y = "-log_10 P-value", title = "BEA Region") + 
  facet_wrap(~region, nrow = 1) + 
  theme_bw() +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 12),  
    axis.title.x = element_text(size = 10),               
    axis.title.y = element_text(size = 10),               
    legend.position = "none"                           
  )




########## state region ################
pval_by_location = readRDS("/Users/sarah/Desktop/ARTP/Real Application/P_individual_state_region.rds")
library(rtp)
library(tidyverse)

### getting ARTP 
ARTP.q = function(p.values) {
  K = length(p.values)
  v.stat.pval = sapply(1:K,function(i) p.rtp(i, p.values))
  v.AFp = min(v.stat.pval)
  return(stat=v.AFp)
}
region_ARTP_P = c()
for (i in 1:length(pval_by_location)) {
  region_p=pval_by_location[i][[1]]
  ARTP_p=c()
  for (j in 1:5) {
    region_p_T = region_p[[j]]
    region_q_T = ARTP.q(region_p_T)
    region_q_T = -log10(region_q_T)
    region_q_T = ifelse(region_q_T>15, 15, region_q_T)
    ### get ARTP p-value
    UFI_base = UFI_database %>% filter(K==length(region_p_T))
    polyd3 <- lm(ARTP_is ~ poly(q, 3), data = UFI_base)
    p_T_UFI <- predict(polyd3, newdata = data.frame(q=region_q_T))
    ARTP_p = c(ARTP_p, p_T_UFI)
  }
  region_ARTP_P = c(region_ARTP_P, list(ARTP_p))
}

## Get Fisher
library(poolr)
region_Fisher_P = c()
for (i in 1:length(pval_by_location)) {
  region_p=pval_by_location[i][[1]]
  Fisher_p=c()
  for (j in 1:5) {
    region_p_T = region_p[[j]]
    Fisher_p_T = (fisher(region_p_T))$p
    Fisher_p = c(Fisher_p, Fisher_p_T)
  }
  region_Fisher_P = c(region_Fisher_P, list(Fisher_p))
}

## Get HC
library(HCp)
region_HC_P = c()
for (i in 1:length(pval_by_location)) {
  region_p=pval_by_location[i][[1]]
  HC_p=c()
  for (j in 1:5) {
    region_p_T = region_p[[j]]
    HC_q_T = HCstat(region_p_T, k0=1, thre=FALSE); HC_q_T
    HC_p_T = ufi.p(flibs=HC_flibs, K=length(pval_by_location), q=HC_q_T); HC_p_T
    HC_p_T <- ifelse(is.nan(HC_p_T), 1, HC_p_T); HC_p_T
    HC_p = c(HC_p, HC_p_T)
  }
  region_HC_P = c(region_HC_P, list(HC_p))
}


names(region_ARTP_P) = names(pval_by_location)
region_plot = as.data.frame(do.call(rbind, region_ARTP_P))
colnames(region_plot) = c("T1","T2","T3","T4","T5"); region_plot$region = rownames(region_plot)

region_plot = pivot_longer(
  region_plot,
  cols = starts_with("T"),        
  names_to = "Period",
  values_to = "-log10P"
)

region_plot$`-log10P` = ifelse(region_plot$`-log10P`<0, 0, region_plot$`-log10P`)
region_plot$ARTP = region_plot$`-log10P`
region_plot$Fisher = ifelse(-log10(unlist(region_Fisher_P))>15, 15, -log10(unlist(region_Fisher_P)))
region_plot$HC = ifelse(-log10(unlist(region_HC_P))>15, 15, -log10(unlist(region_HC_P)))
region_plot = region_plot %>% select(-`-log10P`)

region_plot = pivot_longer(region_plot, cols=c(ARTP, Fisher, HC), 
                           names_to = "method", 
                           values_to = "p_value")

## all midwest state
region_plot %>% filter(region %in% c("Iowa", "Kansas","Minnesota","Missouri","Nebraska","Illinois","Indiana",
                                     "Michigan","Ohio")) %>%
  ggplot(aes(x = Period, y = p_value, color = method,shape = method, group = method)) + 
  geom_point(size=1) + 
  geom_line(alpha=0.5) + 
  labs(x = "Time Period", y = "-log_10 P-value") +
  facet_wrap(~region, nrow = 2) + theme_bw()


## all southeast state
state_region = region_plot %>% filter(region %in% c("Georgia","Kentucky", "Mississippi", "North Carolina",
                                     "Tennessee", "Virginia")) %>% 
  ggplot(aes(x = Period, y = p_value, color = method, shape = method, group = method)) + 
  geom_point(size = 1) + 
  geom_line(alpha = 0.5) + 
  labs(x = "Time Window", y = "-log_10 P-value", title = "States with data from more than 80 counties") + 
  facet_wrap(~region, nrow = 1) + 
  theme_bw() +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 12),  
    axis.title.x = element_text(size = 10),               
    axis.title.y = element_text(size = 10),               
    legend.position = "bottom"                           
  )

## combine all
library(patchwork)
combined_plot <- Geog_region + BEA_region + state_region + 
  plot_layout(ncol = 1)  # stack vertically, change to ncol=3 for side-by-side
combined_plot

ggsave("region_comb.png", plot = combined_plot, width = 7, height = 7) 


