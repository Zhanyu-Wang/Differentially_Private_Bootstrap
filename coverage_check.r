#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

#automatic install of packages if they are not installed already
list.of.packages <- c(
  "foreach",
  "doSNOW",
  "deconvolveR",
  "boot",
  "reshape2"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages) > 0){
  install.packages(new.packages, dep=TRUE)
}

#loading packages
for(package.i in list.of.packages){
  suppressPackageStartupMessages(
    library(
      package.i,
      character.only = TRUE
    )
  )
}

n.cores <- min(124, parallel::detectCores() - 1)
cl <- makeSOCKcluster(n.cores)
registerDoSNOW(cl)

nSIM <- 2000  # simulation number
pb <- txtProgressBar(max = nSIM, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

N <- 10000  # sample size
B_0 <- 1000
pDegree <- 5
c0 <- 1e-1
bin_num <- 1000

xmin <- 0
xmax <- 1
true_mean <- (xmax + xmin)/2
sensitivity <- (xmax - xmin) / N

confidence_level <- 0.9
method_list <- c('deconvolveR')
gdp_mu_list <- c(1, 0.5, 0.3, 0.1)
for (gdp_mu in gdp_mu_list){
  set.seed(12345)
  for (B in c(20, 180, 500, 2000)) {
    print(c(gdp_mu, B))
    all_results <- data.frame(matrix(ncol = 0, nrow = nSIM))
    sd_of_noise <- sqrt(B) * sensitivity / gdp_mu * 1.125
    true_sd_of_sample_mean <- (xmax-xmin)/sqrt(12*N)
    print(c("sd_of_noise, true_sd_of_sample_mean", sd_of_noise, true_sd_of_sample_mean))
    
    
    func_mean <- function(formula, data, indices) {
      d <- data[indices]
      return(mean(d))
    }
    
    func_boot <- function(unif_data) {
      return(boot(data=unif_data, statistic=func_mean, R=B, formula=0)$t) # formula is not used for now
    }
    
    start_time <- Sys.time()
    clean_boot_means <- foreach(
      i = 1:nSIM,
      .combine = 'cbind',
      .packages='boot',
      .options.snow=opts
    ) %dopar% {
      func_boot(pmax(xmin, pmin(xmax, rnorm(n = N, mean=0.5))))
    }
    print(c('bootstrap_step_time', round(difftime(Sys.time(), start_time, unit='secs'), 2)))
    
    
    
    # DP input (noisy results)
    noisy_boot_means <- clean_boot_means + sapply(seq_len(nSIM), function(x) rnorm(n = B, 0, sd_of_noise))
    noisy_boot_means_scaled <- noisy_boot_means / sd_of_noise # from 0.45 to 0.55
    support_min_scaled <- xmin / sd_of_noise
    support_max_scaled <- xmax / sd_of_noise
    bin_width <- (support_max_scaled - support_min_scaled)/bin_num
    grid_seq <- seq(from = support_min_scaled, to = support_max_scaled,
                    by = bin_width)
    bootstrap_step_time <- Sys.time() - start_time
    
    
    # the percentile of true parameter in the recovered distribution
    CDF_to_CI <- function(x_seq, CDF_seq, confidence_level){
      true_param_idx <- which(x_seq > true_mean / sd_of_noise)[1]
      if(is.na(true_param_idx) | true_param_idx==0) true_param_idx=1
      return(CDF_seq[true_param_idx])
    }
    
    func_deconvolveR_CI <- function (x){
      q1 <- quantile(x, 0.25)
      q3 <- quantile(x, 0.75)
      grid_seq <- seq(from=q1-3*(q3-q1), to=q3+3*(q3-q1), length.out=bin_num)
      
      deconv_result <- deconv(tau = grid_seq, X = x, family = "Normal", pDegree = pDegree, c0=c0)
      x_seq <- deconv_result$stats[, 'theta']
      CDF_seq <- deconv_result$stats[, 'G']
      return(CDF_to_CI(x_seq, CDF_seq, confidence_level))
    }
    
    start_time <- Sys.time()
    clean_boot_percentile <- apply(clean_boot_means, 2,
                                   function (clean_t) sum(clean_t < true_mean) / length(clean_t))
    method_name <- paste('clean_bootstrap', sep='')
    all_results[method_name] <- clean_boot_percentile
    
    
    method_name <- paste('deconvolveR', sep='')
    all_results[method_name] <- foreach(
      i = 1:nSIM,
      .combine = 'rbind',
      .packages=c('boot',
                  'deconvolveR'),
      .options.snow=opts
    ) %dopar% {
      func_deconvolveR_CI(noisy_boot_means_scaled[,i])
    }
    
    w.plot <- melt(all_results)
    w.plot <- dplyr::arrange(w.plot, value)
    w.plot$deconvolution_method <- w.plot$variable
    w.plot <- w.plot[order(w.plot['variable']),]
    legend_order <- c(seq(1,ncol(all_results), 2), seq(2,ncol(all_results), 2))
    legend_levels <- c('reference(runif)', colnames(all_results)[legend_order])
    w.plot$deconvolution_method <- factor(w.plot$deconvolution_method,
                                          levels = legend_levels)
    
    csv_name <- paste("results/deconvolution_simulation-", method_name, "_range=", xmax-xmin,
                      "_N=", N, "_mu=", gdp_mu, "_B=", B,
                      "_pd=", pDegree, "_c0=", c0, "_bin=", bin_num,
                      "_nSIM=", nSIM, ".csv", sep='')
    write.csv(w.plot, csv_name, row.names = FALSE)
  }
}
stopCluster(cl)

