
############## Figure Accuracy 1 - showing MC fails for smaller p-values #################
# Run in speed 5
q_input = 10^-(seq(-log10(1e-15),-log10(1), length.out = 15))
K_set = c(10, 50, 100, 200, 500)
idx = ifelse(-log10(q_input)<10, 0.3, 0.5)
N = 1e4

#time_ARTP_1(q=q_input, K=10, method="vse_IS", idx=idx, N=1e2)
vse_MC_p = mclapply(K_set, function(k) {time_ARTP_1(q=q_input, K=k, method="vse_MC", idx=idx, N=N)}, mc.cores=3)
vse_MC_p <- do.call(rbind, vse_MC_p)
vse_IS_p = mclapply(K_set, function(k) {time_ARTP_1(q=q_input, K=k, method="vse_IS", idx=idx, N=N)}, mc.cores=3)
vse_IS_p <- do.call(rbind, vse_IS_p)
vse_UFI_p = mclapply(K_set, function(k) {time_ARTP_1(q=q_input, K=k, method="vse_UFI", idx=idx, N=N)}, mc.cores=3)
vse_UFI_p <- do.call(rbind, vse_UFI_p)

Accuracy_1 = rbind(vse_MC_p, vse_IS_p, vse_UFI_p)
# setwd("/home/yaw130/ARTP")
# saveRDS(Accuracy_1, file="PaperFigure_Accuracy1.rds")
# setwd("/home/yaw130/ARTP")







############ update Accuracy&Timing (Figure5) ###################

# Run in speed 5
q_input = 10^-(seq(-log10(1e-5),-log10(1), length.out = 4))
K_set = c(10, 50, 100, 200, 500)
idx = ifelse(-log10(q_input)<10, 0.3, 0.5)

## start with n=1e3: compare N and computing time that generates 
# the same standard deviation, a concept similar to relative efficiency in estimation.
N = 1e4
vse_MC_p = mclapply(K_set, function(k) {time_ARTP_10(q=q_input, K=k, method="vse_MC", idx=idx, N=N)}, mc.cores=3)
vse_MC_p <- do.call(rbind, vse_MC_p)
vse_IS_p = mclapply(K_set, function(k) {time_ARTP_10(q=q_input, K=k, method="vse_IS", idx=idx, N=N)}, mc.cores=3)
vse_IS_p <- do.call(rbind, vse_IS_p)

vse_MC_p[,c(1, 7:16)] <- lapply(vse_MC_p[,c(1, 7:16)], function(x) ifelse(x == 0, 1 / (1 + N), x))
vse_MC_p[,c(1,7:16)] = -log10(vse_MC_p[,c(1,7:16)])
vse_MC_p$sd = apply(vse_MC_p[, 7:16], 1, function(x) sd(x[x != Inf], na.rm = T))
vse_MC_p$mean = apply(vse_MC_p[, 7:16], 1, function(x) mean(x[x != Inf], na.rm = T))
# vse_MC_p$q = round(vse_MC_p$q, 3)

vse_IS_p[,c(1, 7:16)] <- lapply(vse_IS_p[,c(1, 7:16)], function(x) ifelse(x == 0, 1 / (1 + N), x))
vse_IS_p[,c(1,7:16)] = -log10(vse_IS_p[,c(1,7:16)])
vse_IS_p$sd = apply(vse_IS_p[, 7:16], 1, function(x) sd(x[x != Inf], na.rm = T))
vse_IS_p$mean = apply(vse_IS_p[, 7:16], 1, function(x) mean(x[x != Inf], na.rm = T))
# vse_IS_p$q = round(vse_IS_p$q, 3)


N=1e4
vse_MC_p = mclapply(K_set, function(k) {time_ARTP_10(q=q_input, K=k, method="vse_MC", idx=idx, N=N)}, mc.cores=3)
vse_MC_p <- do.call(rbind, vse_MC_p)

vse_MC_p[,c(1, 7:16)] <- lapply(vse_MC_p[,c(1, 7:16)], function(x) ifelse(x == 0, 1 / (1 + N), x))
vse_MC_p[,c(1,7:16)] = -log10(vse_MC_p[,c(1,7:16)])
vse_MC_p$sd = apply(vse_MC_p[, 7:16], 1, function(x) sd(x[x != Inf], na.rm = T))
vse_MC_p$mean = apply(vse_MC_p[, 7:16], 1, function(x) mean(x[x != Inf], na.rm = T))

Accuracy_Timing = rbind(vse_MC_p, vse_IS_p)
saveRDS(Accuracy_Timing, file="PaperFigure_Accuracy_Timing2.rds")
setwd()

## redo on speed 6



################# Final Plot ##################
method_colors <- c(
  "vse_IS"  = "#1f77b4",  # blue
  "vse_MC"  = "#ff7f0e",  # orange
  "vse_UFI" = "#2ca02c"   # green
)

df = readRDS("result/Timing_Accurancy/PaperFigure_Accuracy1.rds") # K=10,50,100,200,500

p0 <- df %>%
  mutate(p_1 = ifelse(p_1 == 0, 1 / (1 + 1e4), p_1)) %>%
  mutate(q = -log10(q), p_1 = -log10(p_1)) %>%
  ggplot(aes(x = q, y = p_1)) +
  geom_point(aes(color = method), size=0.8,
             position = position_dodge(width = 2)) +
  ylab("-log_10 Estimated P-value") +
  xlab("-log_10 Statistic") +
  theme_bw() +
  facet_wrap(~ K, scales = "free_y", nrow = 1,
             labeller = labeller(K = function(k) paste0("K = ", k))) +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = c(16, 17, 18)) +
  guides(color = guide_legend(title = "Method"))

p0

df = readRDS("result/Timing_Accurancy/PaperFigure_Accuracy_Timing2.rds") %>%
  dplyr::select(K, q, N, method, mean, sd, time)

q_input = 10^-(seq(-log10(1e-5),-log10(1), length.out = 5))
K_set = c(10, 50, 100, 200, 500); N=10

vse_UFI_p = mclapply(K_set, function(k) {time_ARTP_10(q=q_input, K=k, method="vse_UFI", idx=1, N=N)}, mc.cores=1)
vse_UFI_p <- do.call(rbind, vse_UFI_p)
vse_UFI_p[,c(1, 7:16)] <- lapply(vse_UFI_p[,c(1, 7:16)], function(x) ifelse(x == 0, 1 / (1 + N), x))
vse_UFI_p[,c(1,7:16)] = -log10(vse_UFI_p[,c(1,7:16)])
vse_UFI_p$sd = apply(vse_UFI_p[, 7:16], 1, function(x) sd(x[x != Inf], na.rm = T))
vse_UFI_p$mean = apply(vse_UFI_p[, 7:16], 1, function(x) mean(x[x != Inf], na.rm = T))
vse_UFI_p = vse_UFI_p %>% dplyr::select(K, q, N, method, mean, sd, time)
df = rbind(df, vse_UFI_p)


p1 <- df %>%
  ggplot(aes(x = q, y = mean, color = method)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                width = 0.5,
                position = position_dodge(width = 0.8),
                alpha = 0.6) +
  geom_point(position = position_dodge(width = 0.8), size=0.8) +
  ylab("-log_10 Estimated P-value") +
  xlab("-log_10 Statistic") +
  theme_bw() +
  facet_wrap(~ K, scales = "free_y", nrow = 1,
             labeller = labeller(K = function(k) paste0("K = ", k))) +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = c(16, 17, 18)) +
  guides(color = guide_legend(title = "Method"))
p1


p2 <- df %>%
  filter(method != "vse_MC") %>%
  ggplot(aes(x = K, y = time, color = method, group = method)) +
  geom_line() +
  geom_point(size=0.8) +
  ylab("Time (second)") +
  xlab("K") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  facet_wrap(~ q, scales = "free_y", nrow = 1,
             labeller = labeller(q = function(q) paste0("-log_10 (r) = ", q))) +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = c(16, 17, 18)) +
  guides(color = guide_legend(title = "Method"))


library(ggplot2)
library(patchwork)
p = p0 / p1 / p2  
p
saveRDS(p, file="PaperFigure_Accuracy_Timing.rds")






############ Figure Timing - comparing timing for MC, ISCE, UFI ####################
q.val = 10^-(seq(-log10(1e-15),-log10(1), length.out = 15))
q_input = q.val
K_set = round(seq(5,500, length.out = 20))
idx = ifelse(-log10(q_input)<10, 0.3, 0.5)
N = 1e4

#time_ARTP_1(q=q_input, K=10, method="vse_IS", idx=idx, N=1e2)
vse_MC_p = mclapply(K_set, function(k) {time_ARTP_1(q=q_input, K=k, method="vse_MC", idx=idx, N=N)}, mc.cores=3)
vse_MC_p <- do.call(rbind, vse_MC_p)
vse_IS_p = mclapply(K_set, function(k) {time_ARTP_1(q=q_input, K=k, method="vse_IS", idx=idx, N=N)}, mc.cores=3)
vse_IS_p <- do.call(rbind, vse_IS_p)
vse_UFI_p = mclapply(K_set, function(k) {time_ARTP_1(q=q_input, K=k, method="vse_UFI", idx=idx, N=N)}, mc.cores=3)
vse_UFI_p <- do.call(rbind, vse_UFI_p)

timing_1e4 = rbind(vse_MC_p, vse_IS_p, vse_UFI_p)
setwd("/home/yaw130/ARTP")
save(timing_1e5, file="timing_1e5.RData")
setwd("/home/yaw130/ARTP")


















