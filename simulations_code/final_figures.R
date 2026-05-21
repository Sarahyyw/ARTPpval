
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)


##### Figure 1 -  ISCE Performance with five different proposal #####
load("result/ARTP/previous/ARTP_50_all_dist.RData")

p = res_all_dist %>% filter(q>3) %>%
  filter(idx %in% c("X", 1), N == 1e4) %>% 
  ggplot(aes(x = q, y = mean, group = family, color = family)) + 
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, color = family, width = 0.5), 
                position = position_dodge(width = 0.8), width=0.3, size=1) + 
  geom_point(aes(color = family, shape = family), 
             position = position_dodge(width = 0.8), size=3) + 
  ylab("-log_10 Estimated P-value") + 
  xlab("log_10 statistic") + 
  theme_bw() +  # White background theme
  # ggtitle("ISCE Performance with five different proposal \n density families for K = 50") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),  
    axis.text.y = element_text(size = 14),  
    axis.title.x = element_text(size = 16, face = "bold"),  
    axis.title.y = element_text(size = 16, face = "bold"), 
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
    legend.text = element_text(size = 14),  
    legend.title = element_text(size = 16), 
    legend.position = "bottom"
  )
# ggsave("plot.png", plot = p, width = 10, height = 6) 





##### Figure 2 - UFI interpolation #####
start <- 1; end <- 1000; increment <- 100
ARTP_res <- c()
for (i in seq(start, end, by = increment)) {
  file_name <- paste0("result/ARTP/UFI/ARTP_", i, "_", i + increment - 1, "_res_grad_1e4.RData")
  load(file_name)
  ARTP_res <- c(ARTP_res, get(paste0("ARTP_", i, "_", i + increment - 1, "_res_grad_1e4")))
}

combined_list <- lapply(ARTP_res, function(x) {return(x$res)})
res <- do.call(rbind, combined_list)
res[,c(1,7)] = -log10(res[,c(1,7)])

load("result/ARTP/UFI/ARTP_test.RData")
test_data[,c(2,5)] = -log10(test_data[,c(2,5)])
test_data$pred_d1=NA; test_data$pred_d2=NA; test_data$pred_d3=NA; test_data$pred_d4=NA

for (r in 1:nrow(test_data)) {
  k=test_data$K[r]; q=test_data$q[r]
  train = res %>% filter(K==k)
  polyd1 <- lm(ARTP_is ~ poly(q, 1), data = train)
  test_data$pred_d1[r] <- predict(polyd1, newdata = data.frame(q=q))
  
  polyd2 = lm(ARTP_is ~ poly(q, 2), data = train)
  test_data$pred_d2[r] <- predict(polyd2, newdata = data.frame(q=q))
  
  polyd3 = lm(ARTP_is ~ poly(q, 3), data = train)
  test_data$pred_d3[r] <- predict(polyd3, newdata = data.frame(q=q))
  
  polyd4 = lm(ARTP_is ~ poly(q, 4), data = train)
  test_data$pred_d4[r] <- predict(polyd4, newdata = data.frame(q=q))
}

test_data$p = test_data$ARTP_is
test_data$d1_diff = test_data$pred_d1-test_data$p
test_data$d2_diff = test_data$pred_d2-test_data$p
test_data$d3_diff = test_data$pred_d3-test_data$p
test_data$d4_diff = test_data$pred_d4-test_data$p

# Create the individual violin plots as before
common_theme <- theme_bw() +
  theme(
    axis.text.x  = element_text(size = 14, face = "bold"),
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )

d1 <- test_data %>%
  ggplot() +
  geom_violin(aes(x = "linear", y = d1_diff), color = "blue", fill = "blue") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylim(-0.4, 0.4) +
  labs(x = "") +
  common_theme

d2 <- test_data %>%
  ggplot() +
  geom_violin(aes(x = "quadratic", y = d2_diff), color = "red", fill = "red") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylim(-0.4, 0.4) +
  labs(x = "") +
  common_theme

d3 <- test_data %>%
  ggplot() +
  geom_violin(aes(x = "Cubic", y = d3_diff), color = "green", fill = "green") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylim(-0.4, 0.4) +
  labs(x = "") +
  common_theme

d4 <- test_data %>%
  ggplot() +
  geom_violin(aes(x = "Bi-quadratic", y = d4_diff), color = "blue", fill = "blue") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylim(-0.4, 0.4) +
  labs(x = "") +
  common_theme

library(cowplot)
p_main <- plot_grid(d1, d2, d3, d4, nrow = 1)
y_lab <- ggdraw() +
  draw_label("Log relative error (interpolation vs ISCE)",
             angle = 90, size = 14)
p <- plot_grid(y_lab, p_main, nrow = 1, rel_widths = c(0.08, 1))
p
# ggsave("plot.png", plot = p, width = 6, height = 3) 








##### Figure 3 - Prob bound #####
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

K_vals <- c(5, 10, 50, 100, 200, 300)   # <-- change these to your five K values
beta_seq <- seq(from = 0, to = 1, by = 0.1)
j_levels <- c("S", "2S")

make_df_for_K <- function(K) {
  rtp_K <- readRDS(paste0("result/RTP/RTP_noLog_tau_K=", K, ".rds"))
  # rtp_K <- readRDS(paste0("result/RTP/withoutLog/RTP_noLog_K=", K, ".rds"))
  
  # keep columns: num_sig (1), gamma (2), check_point (3), N (21)
  res_heatmap <- cbind(rtp_K$res, rtp_K$res_bench)[, c(1, 2, 3, 6, 19, 20)]
  
  # assign beta and j (ensure lengths match rows)
  res_heatmap$bound_fail = res_heatmap$benchmark > sapply(res_heatmap$p_rtp, benchmark_tau)
  res_heatmap$bound = sapply(res_heatmap$p_rtp, benchmark_tau)
  res_heatmap$p_rtp = res_heatmap$p_rtp
  res_heatmap$beta <- rep(beta_seq, each = 6)  # matches your earlier pattern
  res_heatmap$j    <- rep(j_levels, length.out = nrow(res_heatmap))
  
  res_heatmap = res_heatmap %>%
    mutate(
      K        = factor(paste0("K = ", K), levels = paste0("K = ", K_vals)),
      beta     = factor(beta),
      gamma    = factor(gamma),
      check_point = factor(j, levels = j_levels),
      gamma_j  = paste0("γ=", gamma, ", j=", check_point),
      N        = factor(N, levels = c("10000", "1e+05", "1e+06")),
      # bound_success = as.logical(bound_success)
    ) %>%
    select(K, beta, gamma_j, N, bound_fail, p_rtp, benchmark, bound)
}

df_all <- map_dfr(K_vals, make_df_for_K)


y_labels <- c(
  "γ=1e-07, j=S"  = expression(gamma == 10^{-7} * ", j = S"),
  "γ=1e-07, j=2S" = expression(gamma == 10^{-7} * ", j = 2S"),
  "γ=1e-05, j=S"  = expression(gamma == 10^{-5} * ", j = S"),
  "γ=1e-05, j=2S" = expression(gamma == 10^{-5} * ", j = 2S"),
  "γ=0.001, j=S"  = expression(gamma == 10^{-3} * ", j = S"),
  "γ=0.001, j=2S" = expression(gamma == 10^{-3} * ", j = 2S"))

p <- ggplot(df_all, aes(x = beta, y = gamma_j, fill = N)) +
  geom_tile(color = "white") +
  scale_fill_manual(
    name = expression(N),
    breaks = c("10000", "1e+05", "1e+06"),
    labels = c("10000" = expression(10^4),
               "1e+05" = expression(10^5),
               "1e+06" = expression(10^6)),
    values = c("10000" = "lightblue", "1e+05" = "yellow", "1e+06" = "orange"),
    na.value = "grey90") +
  scale_y_discrete(labels = y_labels) +
  facet_wrap(~ K, nrow = 2) +
  labs(
    # title = "Heatmaps for benchmarking RTP",
    x = expression(beta),
    y = expression(gamma * ", j (combined)")
  ) +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=7),
        axis.text.x = element_text(angle=70, hjust=0.8, size=7),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold")
  )

p 
# ggsave("RTP_benchmark.png", plot = p, width = 6, height = 3) 





##### Figure 4 - time & accuracy #####
## run the Timming_server.R first
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
#ggsave("PaperFigure_Accuracy_Timing..png", plot = p, width = 6, height = 6) 





##### Figure 5 - Application #####


















