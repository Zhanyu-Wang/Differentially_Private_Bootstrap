#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# ind_data <- read.csv(file = 'ind2016.csv')
# ind_data2 <- ind_data[,c(94,115,130)]
# write.csv(ind_data2, 'ind2016_MRKINC-PR-SHELCO.csv')
ind_data2 <- read.csv(file = 'ind2016_MRKINC-PR-SHELCO.csv')
nSIM <- 2000
result_digits <- 10

# provide_PR_code_name <- "Quebec"
provide_PR_code_name <- "Ontario"

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

description_str <- "ind2016-logistic_"

row_col <- dim(ind_data2) # 930421    4

# col_index:  MRKINC=94, SHELCO=130
# PR:         Quebec=24, Ontario=35
provide_PR_codes <- data.frame(Quebec=24, Ontario=35)
col_indices <- data.frame(MRKINC=94, SHELCO=130)        

provide_PR_code <- as.numeric(provide_PR_codes[provide_PR_code_name])
province_data <- ind_data2[ind_data2[1:row_col[1], "PR"] == provide_PR_code, colnames(col_indices)]
row_indice <- (province_data != 99999999) & (province_data != 88888888) &
  (province_data != 678800) & (province_data %% 100 == 0) & (province_data >= 0)
row_indice_1 <- row_indice[,1] & row_indice[,2]
province_data <- province_data[row_indice_1,]
# print(c(max(province_data[,1]), max(province_data[,2])))
max_x <- max(province_data[,1])
max_y <- max(province_data[,2])
province_data[,1] <- province_data[,1] / max(province_data[,1])
province_data[,2] <- province_data[,2] / max(province_data[,2])
population_size <- dim(province_data)[1]


expit = function(x){
  return(exp(x)/(1+exp(x)))
}

LogisticRegularized = function(X, Y, penalty){
  n = length(Y)
  obj = function(beta){
    return((1/n)*sum(log(1+exp(-Y * X%*%beta))) + penalty*sum(beta^2))
  }
  grad = function(beta){
    return( -(1/n) * t(Y * expit(-Y * X%*%beta)) %*% X + 2*penalty*beta)
  }
  min = optim(par = rep(0, 2), fn = obj, gr = grad, method = "BFGS")
  
  if(min$convergence!=0){
    print("Output-Pert did not converge") 
  }
  return(min$par)
}###   END OBJ  PERT



PrivSPDMat = function(original_matrix, sensitivity, ep, c){
  d = dim(original_matrix)[1]
  vector = rnorm(n=d^2, m=0, s=1)
  eta = matrix(vector * sensitivity / ep, nrow=d)
  new_matrix = original_matrix + eta
  new_matrix = (new_matrix + t(new_matrix)) / 2
  ev <- eigen(new_matrix)
  new_ev <- pmax(2*c, ev$values)
  new_matrix2 <- ev$vectors %*% diag(new_ev) %*% t(ev$vectors)
  return(new_matrix2)
}


###################################
   
confidence_levels <- c(0.9, 0.95)

for (penalty_c in c(0.01, 1)){ 
  for (gdp_mu in c(1, 2, 5, 10)){
    for (B in c(20,50,100,200,500,1000)){
      for (N in c(1000, 3000, 10000, 30000, 100000)){ # 
        print(c(N, penalty_c, gdp_mu))
        txt_name <- paste("results_logistic/deconvolution_", description_str, provide_PR_code_name, 
                          "_N=", N, "_mu=", gdp_mu, "_B=", B, "_c=", penalty_c, "_nSIM=", nSIM, ".csv", sep='')
        
        if (file.exists(txt_name)) {
          next
        }
        X <- cbind(1, province_data$MRKINC) / sqrt(2)
        Y <- as.integer(province_data$SHELCO >= 0.5) * 2 - 1

        true_param <- LogisticRegularized(X, Y, penalty=penalty_c)[2]
        print(true_param)
        ############################################################
        # This experiment focuses on logistic
        func_param <- function(formula, data, indices) {
          X <- cbind(1, data[indices, 1]) / sqrt(2)
          Y <- as.integer(data[indices, 2] >= 0.5) * 2 - 1
          return(LogisticRegularized(X, Y, penalty=penalty_c)[2])
        }

        sensitivity <- 1 / N / penalty_c
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

        func_boot <- function(unif_data) {
          dim(unif_data) <- c(N,2)
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
          func_boot(as.matrix(province_data[sample(1:population_size, N, replace=TRUE),]))
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
        method_name <- 'DP-CI-ERM'
        start_time <- Sys.time()

        method_result <- foreach(
          i = 1:nSIM, 
          .combine = 'cbind',
          .packages='boot',
          .options.snow=opts
        ) %dopar% {
          set.seed(i)
          unif_data <- province_data[sample(1:population_size, N, replace=TRUE), ] 

          ep1 = gdp_mu/2
          ep2 = ep1/2
          ep3 = ep1/2

          l2sensitivity = 1 / N / penalty_c
          noise_level = l2sensitivity / ep1
          
          X <- cbind(1, unif_data[, 1]) / sqrt(2)
          Y <- as.integer(unif_data[, 2] >= 0.5) * 2 - 1
              
          x_dim <- dim(X)[2]-1
          nonprivate_beta = LogisticRegularized(penalty=penalty_c, X=X, Y=Y)
          private_beta = nonprivate_beta + rnorm(n=1, m=0, s=1) * noise_level
          
          new_matrix = X
          multiplier = expit(Y * X%*%private_beta) * expit(-Y * X%*%private_beta)
          new_matrix[,1] = multiplier * new_matrix[,1]
          new_matrix[,2] = multiplier * new_matrix[,2]
          hessian_estimate = t(new_matrix) %*% X / N + 2 * penalty_c * diag(2)
          print("hessian_estimate")
          print(hessian_estimate)
          hessian_sensitivity = 1/(2*N)
          private_hessian = PrivSPDMat(hessian_estimate, hessian_sensitivity, ep2, penalty_c)
          private_hessian_inv = solve(private_hessian)
          print("private_hessian")
          print(private_hessian)
          print(private_hessian_inv)
          
          new_matrix = X
          multiplier = expit(-Y * X%*%private_beta)^2
          new_matrix[,1] = multiplier * new_matrix[,1]
          new_matrix[,2] = multiplier * new_matrix[,2]
          cov_estimate = t(new_matrix) %*% X / N - 4 * penalty_c^2 * private_beta %*% t(private_beta)
          print("cov_estimate")
          print(cov_estimate)
          cov_sensitivity = 2/N * expit(sqrt(sum(private_beta^2)))^2
          ### the original algorithm in DP-CI-ERM used penalty_c for private covariance estimation, but it is not necessary.
          # private_cov = PrivSPDMat(cov_estimate, cov_sensitivity, ep3, penalty_c)   
          private_cov = PrivSPDMat(cov_estimate, cov_sensitivity, ep3, 0) 
          print("private_cov")
          print(private_cov)

          U = diag(2) * noise_level^2 + 1 / N * private_hessian_inv %*% private_cov %*% private_hessian_inv
          print(diag(2) * noise_level^2)
          print("U")
          print(U)
            
          method_result_new <- c()
          for (confidence_level in confidence_levels) {
            method_result_new <- rbind(method_result_new, 
                                      private_beta[2] - qnorm((1+confidence_level)/2) * sqrt(U[2,2]),
                                      private_beta[2] + qnorm((1+confidence_level)/2) * sqrt(U[2,2]))
          }
          print("nonprivate_beta")
          print(nonprivate_beta)
          print(private_beta)
          print(method_result_new)
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
        print(txt_name)
        write.table(all_results, txt_name, sep=",", row.names=FALSE, col.names=TRUE)
      }
    }
  }
}