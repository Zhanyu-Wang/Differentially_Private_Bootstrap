#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

nSIM <- 2000
result_digits <- 10

#automatic install of packages if they are not installed already
list.of.packages <- c(
  "foreach",
  "doSNOW",
  "deconvolveR",
  "MASS",
  "boot"
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

pb <- txtProgressBar(max = nSIM, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)



###################################
for (confidence_level in c(0.9)){
  confidence_levels <- c(confidence_level)
  for (gdp_mu in c(0.1, 0.2, 0.5, 1)){
    for (N in c(100, 300, 1000, 3000, 10000, 30000, 100000)){ # 
      for (B in c(1000, 500, 200, 100, 50, 30, 20, 10, 5)){
  
        true_param <- 0.5
        sensitivity <- 1 / N
        sd_of_noise <- sqrt(B) * sensitivity / gdp_mu * 1.125
  
        ############################################################
  
        # for each experiment, we care about the coverage and the width
        CDF_to_CI <- function(x_seq, CDF_seq, confidence_levels){
          true_param_idx <- which(x_seq > true_param / sd_of_noise)[1]
          if(is.na(true_param_idx) | true_param_idx==0) true_param_idx=1
          percentile_of_true_val <- CDF_seq[true_param_idx]
  
          true_param_idx2 <- which(x_seq > 0 / sd_of_noise)[1]
          if(is.na(true_param_idx2) | true_param_idx2==0) true_param_idx2=1
          percentile_of_true_val2 <- CDF_seq[true_param_idx2]
          # this percentile indicates coverage
          output_list <- c(percentile_of_true_val, percentile_of_true_val2)
          
          for (confidence_level in confidence_levels) {
            lowerbound_idx <- which(CDF_seq > (1-confidence_level)/2)[1]-1
            upperbound_idx <- which(CDF_seq > (1+confidence_level)/2)[1]
            if(lowerbound_idx==0) lowerbound_idx=1
            if(upperbound_idx==0) upperbound_idx=length(CDF_seq)
            lowerbound <- x_seq[lowerbound_idx]
            upperbound <- x_seq[upperbound_idx]
            output_list <- c(output_list, upperbound-lowerbound)
          }
          return(output_list)
        }
        bootestimates_to_CI <- function(bootestimates, confidence_levels){
          percentile_of_true_val <- sum(bootestimates <= true_param) / length(bootestimates)
          percentile_of_true_val2 <- sum(bootestimates <= 0) / length(bootestimates)
          # this percentile indicates coverage
          output_list <- c(percentile_of_true_val, percentile_of_true_val2)
          
          for (confidence_level in confidence_levels) {
            lowerbound <- quantile(bootestimates, (1-confidence_level)/2)
            upperbound <- quantile(bootestimates, (1+confidence_level)/2)
            output_list <- c(output_list, upperbound-lowerbound)
          }
          return(output_list)
        }
  
        CI_coverage_width <- function(CI, method_name=NA){
          CI_infos <- c()
          conf_lvl_index <- 2
          for (confidence_level in confidence_levels) {
            conf_lvl_index <- conf_lvl_index + 1
            coverage <- (CI[1,] <= (1+confidence_level)/2) & (CI[1,] >= (1-confidence_level)/2)
            coverage2 <- (CI[2,] <= (1+confidence_level)/2) & (CI[2,] >= (1-confidence_level)/2)
            mean_coverage <- mean(coverage)
            mean_coverage2 <- mean(coverage2)
            # std_coverage <- round(sd(coverage) / sqrt(nSIM), digits = result_digits)
            width <- CI[conf_lvl_index,]
            mean_width <- round(mean(width), digits = result_digits)
            std_width <- round(sd(width) / sqrt(nSIM), digits = result_digits)
            # CI_infos <- c(CI_infos, mean_coverage, std_coverage, mean_width, std_width)
            CI_infos <- c(CI_infos, mean_coverage, mean_coverage2, mean_width, std_width)
          }
          result <- data.frame(method_name, data.frame(matrix(CI_infos, nrow=1)))
          return(result)
        }
  
  
        ### deconvolveR Bootstrap DP CI
        pDegree <- 5
        c0 <- 1e-1
        bin_num <- 1000
        func_deconvolveR_CI <- function (x){
          q1 <- quantile(x, 0.25)
          q3 <- quantile(x, 0.75)
          # make grid for the range of calculating the deconvoluted results
          grid_seq <- seq(from=q1-3*(q3-q1), to=q3+3*(q3-q1), length.out=bin_num) 
          deconv_result <- deconv(tau = grid_seq, X = x, family = "Normal", pDegree = pDegree, c0=c0)
          x_seq <- deconv_result$stats[, 'theta']
          CDF_seq <- deconv_result$stats[, 'G']
          return(CDF_to_CI(x_seq, CDF_seq, confidence_levels))
        }
    
        func_param <- function(formula, data, indices) {
          d <- data[indices]
          return(mean(d))
        }
  
        func_boot <- function(unif_data) {
          # dim(unif_data) <- c(N,2)
          return(boot(data=unif_data, statistic=func_param, R=B, formula=0)$t) # formula is not used for now
        }
  
        set.seed(42)
  
        start_time <- Sys.time()
        clean_boot_params <- foreach(
          i = 1:nSIM, 
          .combine = 'cbind',
          .packages='boot',
          .options.snow=opts
        ) %dopar% {
          set.seed(i)
          func_boot(pmax(0, pmin(1, rnorm(n = N, mean=0.5))))
        }
        bootstrap_step_time <- Sys.time() - start_time
  
  
        ####################################
        method_name <- 'clean_bootstrap'
        start_time <- Sys.time()
  
        clean_boot_CI <- foreach(
          i = 1:nSIM, 
          .combine = 'cbind',
          .packages='boot',
          .options.snow=opts
        ) %dopar% {
          bootestimates_to_CI(clean_boot_params[,i], confidence_levels)
        }
  
        method_result <- clean_boot_CI
        method_stat <- CI_coverage_width(method_result, method_name)
        method_stat$time <- round(Sys.time() - start_time + bootstrap_step_time, 2)
        all_results <- method_stat
  
        ####################################
        method_name <- 'deconvolveR'
        start_time <- Sys.time()
  
        # DP input (noisy results)
        noisy_boot_params <- clean_boot_params + sapply(seq_len(nSIM), function(x) rnorm(n = B, 0, sd_of_noise))
        noisy_boot_params_scaled <- noisy_boot_params / sd_of_noise # from 0.45 to 0.55
        support_min_scaled <- 0 / sd_of_noise 
        support_max_scaled <- 1 / sd_of_noise 
        bin_width <- (support_max_scaled - support_min_scaled)/bin_num
        grid_seq <- seq(from = support_min_scaled, to = support_max_scaled, 
                        by = bin_width)
  
  
        method_result <- foreach(
          i = 1:nSIM, 
          .combine = 'cbind',
          .packages=c('boot',
                      'deconvolveR'),
          .options.snow=opts
        ) %dopar% {
          set.seed(i)
          func_deconvolveR_CI(noisy_boot_params_scaled[,i])
        }
  
        method_result[3:(2+length(confidence_levels)),] <- 
          method_result[3:(2+length(confidence_levels)),] * sd_of_noise
        method_stat <- CI_coverage_width(method_result, method_name)
        method_stat$time <- round(Sys.time() - start_time + bootstrap_step_time, 2)
        all_results <- rbind(all_results, method_stat)
  
  
        ####################################
  
  
        ####################################
        method_name <- 'NoisyVar'
        start_time <- Sys.time()
  
        clean_means_vars <- foreach(
          i = 1:nSIM, 
          .combine = 'cbind',
          .packages='boot',
          .options.snow=opts
        ) %dopar% {
          set.seed(i)
          curr_data <- (pmax(0, pmin(1, rnorm(n = N, mean=0.5))))
          result <- rbind(mean(curr_data), var(curr_data))
        }
        clean_means <- clean_means_vars[1,]
        clean_vars <- clean_means_vars[2,]
  
        sd_of_noise_mean <- sqrt(2) * sensitivity / gdp_mu
        noisy_means <- clean_means + rnorm(n = nSIM, 0, sd_of_noise_mean)
        sensitivity_var <- 1 / N
        sd_of_noise_var <- sqrt(2) * sensitivity_var / gdp_mu
        noisy_vars <- clean_vars + rnorm(n = nSIM, 0, sd_of_noise_var)
        noisy_sds <- sqrt(pmax(noisy_vars / N, 0.0000000001)) 
  
        ##### SIM #####
        nsim_sub <- 1000
        method_result <- foreach(
          i = 1:nSIM, 
          .combine = 'cbind',
          .options.snow=opts
        ) %dopar% {
          clean_means_new <- sapply(seq_len(nsim_sub), function(x) mean(pmin(1, pmax(0, rnorm(n = N, mean=noisy_means[i], sd=noisy_sds[i]*sqrt(N))))))
          noisy_means_new <- clean_means_new + rnorm(n = nsim_sub, 0, sd_of_noise_mean)
          method_result_new <- c()
          for (confidence_level in confidence_levels) {
            moe <- quantile(noisy_means_new, probs=c((1-confidence_level)/2, (1+confidence_level)/2), names=FALSE)
            method_result_new <- rbind(method_result_new, noisy_means[i] - (moe[2]-moe[1])/2,
                                      noisy_means[i] + (moe[2]-moe[1])/2)
          }
          method_result_new
        }
  
        CI_infos <- c()
        for (i in c(1:length(confidence_levels))) {
          coverage <- (method_result[2*i-1,] <= true_param) & (method_result[2*i,] >= true_param)
          mean_coverage <- mean(coverage)
          coverage2 <- (method_result[2*i-1,] <= 0) & (method_result[2*i,] >= 0)
          mean_coverage2 <- mean(coverage2)
          width <- (c(method_result[2*i,] - method_result[2*i-1,])) 
          mean_width <- round(mean(width), digits = result_digits)
          std_width <- round(sd(width) / sqrt(nSIM), digits = result_digits)
          CI_infos <- c(CI_infos, mean_coverage, mean_coverage2, mean_width, std_width)
        }
        method_stat <- data.frame(method_name, data.frame(matrix(CI_infos, nrow=1)))
        method_stat$time <- round(Sys.time() - start_time, 2)
        print(method_stat)
        all_results <- rbind(all_results, method_stat)
  
  
        ####################################
        ##### finish running experiment; print and save results
        print(all_results)
        colname_all_results <- c("method")
        for (confidence_level in confidence_levels) {
          colname_all_results <- c(colname_all_results, 
                                  paste("CI coverage (mean): level=", confidence_level, sep=''),
                                  paste("CI coverage2 (mean): level=", confidence_level, sep=''),
                                  paste("CI width (mean): level=", confidence_level, sep=''),
                                  paste("CI width (se): level=", confidence_level, sep=''))
        }
        colnames(all_results) <- c(colname_all_results, "time")
        txt_name <- paste("results_clamp_normal_mean/deconvolution_clamp_normal_mean_N(0.5,1)_[0,1]_", 
                          "_N=", N, "_mu=", gdp_mu, "_B=", B, "_nSIM=", nSIM, "_", confidence_level, ".csv", sep='')
        print(txt_name)
        write.table(all_results, txt_name, sep=",", row.names=FALSE, col.names=TRUE)
      }
    }
  }
}