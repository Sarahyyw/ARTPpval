

# -----------------  START HERE  -------------------
## K=5 (speed 5 done)
## K=10 (speed 5 done)
## K=50 (speed 6 done)
## K=100 (speed 8 done)
## K=200 (speed 5 done)
## K=300 (speed 9 done)
## K=500 (drop this)

K = 300
# set.seed(697069) 
beta =  seq(from = 0, to = 1, by = 0.1)
gamma = c(10^(-3), 10^(-5), 10^(-7)) ## target p-value

p_input_K = input_generate(K, beta, gamma)
p_input = p_input_K
sig_sum = p_input$signal_sum; sig_unique = unique(sig_sum$num_sig)

res = run_for_RTP(p_input = p_input, N = 1e4, n_input=K)
res_bench = benchmark_func(res) 

## check how many we need to re-run
res_bench$benchmark <= sapply(res_bench$p_rtp, benchmark_tau)



# -----------------  redo the row by increase N (without log) -------------------
for (i in seq_len(nrow(res_bench))) {
  benchmark_i <- res_bench$benchmark[i]
  p_rtp_i <- res_bench$p_rtp[i]  
  bound_i  <- benchmark_tau(p_rtp_i)
  
  if (benchmark_i > bound_i) {
    N <- 1e5
    res_redo_row <- process_row(
      i, res = res, N = N, p_input = p_input,
      n_input = K, sig_unique = sig_unique, gamma = gamma
    ) ## im here
    res_bench_redo_row <- benchmark_func(res_redo_row)
    benchmark_i <- res_bench_redo_row$benchmark
    
    # If still above bound, try N = 1e6 with a hard cap on attempts
    attempts_1e6 <- 0L
    max_attempts_1e6 <- 3L  # hard stop 
    
    while (benchmark_i > bound_i && attempts_1e6 < max_attempts_1e6) {
      N <- 1e6
      
      # After 3 attempts at 1e6, regenerate p_input each try
      if (attempts_1e6 >= 1) {p_input <- input_generate(K, beta, gamma)}
      
      res_redo_row <- process_row(
        i, res = res, N = N, p_input = p_input,
        n_input = K, sig_unique = sig_unique, gamma = gamma
      )
      res_bench_redo_row <- benchmark_func(res_redo_row)
      benchmark_i <- res_bench_redo_row$benchmark
      attempts_1e6 <- attempts_1e6 + 1L
    }
    
    if (benchmark_i > bound_i) {
      warning(sprintf(
        "Row %d: still above bound after %d attempt(s) at N=1e6; keeping last result.",
        i, attempts_1e6))
    }
    
    message(sprintf(
      "redo_done_row = %d | N = %s | benchmark = %.6f | bound = %.6f | attempts_1e6 = %d",
      i, format(N, scientific = FALSE), round(benchmark_i, 6), round(bound_i, 6), attempts_1e6
    ))
    
    res[i, ]       <- res_redo_row
    res_bench[i, ] <- res_bench_redo_row
  }
}

## K=300 (speed5 and speed9)
res_final=list(res=res, res_bench=res_bench)
# saveRDS(res_final, file = "RTP_noLog_tau_K=300.rds")



# -----------------  Res Heatmap  -------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

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

K_vals <- c(5, 10, 50, 100, 200,300)   # <-- change these to your five K values
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
  "γ=1e-07, j=2S" = expression(gamma == 10^{-7} * ", j = 2*S"),
  "γ=1e-05, j=S"  = expression(gamma == 10^{-5} * ", j = S"),
  "γ=1e-05, j=2S" = expression(gamma == 10^{-5} * ", j = 2*S"),
  "γ=0.001, j=S"  = expression(gamma == 10^{-3} * ", j = S"),
  "γ=0.001, j=2S" = expression(gamma == 10^{-3} * ", j = 2*S"))

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
    title = "Heatmaps for benchmarking RTP",
    x = expression(beta),
    y = expression(gamma * ", j (combined)")
  ) +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=8),
        axis.text.x = element_text(angle=70, hjust=0.8, size=8),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold")
        )
p



# rtp_K <- readRDS(paste0("result/RTP/withoutLog/RTP_noLog_K=", 200, ".rds"))
# res = rtp_K$res; res_bench = rtp_K$res_bench
# 
# a <- cbind(rtp_K$res, rtp_K$res_bench)
# a <- a[, !duplicated(colnames(a))]
# a <- dplyr::select(a, num_sig, p_rtp, N, benchmark)
# b = a %>% filter(benchmark > 0.1)
# 



final_res = readRDS("RTP_noLog_tau_K=200.rds")
res = final_res$res
res_bench = final_res$res_bench



res_bench$benchmark <= sapply(res_bench$p_rtp, benchmark_tau)
FALSE_row = which(!res_bench$benchmark <= sapply(res_bench$p_rtp, benchmark_tau)); FALSE_row 

for (i in FALSE_row) {
  N = 1e6
  p_input <- input_generate(K, beta, gamma)
  res_redo_row <- process_row(
    i, res = res, N = N, p_input = p_input,
    n_input = K, sig_unique = sig_unique, gamma = gamma
  )
  benchmark_func(res_redo_row)
  
  res_bench_redo_row <- benchmark_func(res_redo_row)
  benchmark_i <- res_bench_redo_row$benchmark
  bound_i  <- benchmark_tau(res_bench_redo_row$p_rtp)
  message(sprintf(
    "redo_done_row = %d | N = %s | benchmark = %.6f | bound = %.6f | attempts_1e6 = %d",
    i, format(N, scientific = FALSE), round(benchmark_i, 6), round(bound_i, 6), attempts_1e6
  ))
  
  res[i, ]       <- res_redo_row
  res_bench[i, ] <- res_bench_redo_row
}





