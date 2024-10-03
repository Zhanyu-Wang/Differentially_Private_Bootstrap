CDF_to_CI <- function(x_seq, CDF_seq, confidence_level){
  true_param_idx <- which(x_seq > true_mean / sd_of_noise)[1]
  if(is.na(true_param_idx) | true_param_idx==0) true_param_idx=1
  return(CDF_seq[true_param_idx])
}

library(boot)
library(ggplot2)
library(reshape2)

N <- 10000
B <- 1000
gdp_mu <- 1

nSIM <- 10  # simulation number

xmin <- 0
xmax <- 1
confidence_level <- 0.9
bin_num <- 300
all_results <- data.frame(matrix(ncol = 0, nrow = nSIM))
set.seed(12346)

true_mean <- (xmax + xmin)/2
sensitivity <- (xmax - xmin) / N
sd_of_noise <- sqrt(B) * sensitivity / gdp_mu
nSIM_data <- sapply(seq_len(nSIM), function(x) runif(n = N, min=xmin, max=xmax))
true_sd_of_sample_mean <- (xmax-xmin)/sqrt(12*N)

print(c("sd_of_noise, true_sd_of_sample_mean", sd_of_noise, true_sd_of_sample_mean))

func_mean <- function(formula, data, indices) {
d <- data[indices] 
return(mean(d))
}

func_boot <- function(unif_data) {
return(boot(data=unif_data, statistic=func_mean, R=B, formula=0)$t) # formula is not used for now
}

clean_boot_means <- apply(nSIM_data, 2, func_boot)
# DP input (noisy results)
start_time <- Sys.time()
noisy_boot_means <- clean_boot_means + sapply(seq_len(nSIM), function(x) rnorm(n = B, 0, sd_of_noise))
noisy_boot_means_scaled <- noisy_boot_means / sd_of_noise # from 0.45 to 0.55
bootstrap_step_time <- Sys.time() - start_time


library(deconvolveR) # only support normal

den_x <- c()
den_y <- c()
facet_idx <- c()
empirical_density_color <- c()
empirical_density_linetype <- c()
for (i in seq(3)){
  c0 <- 1e-1
  test_method <- 'deconvolveR'
  
  curr_data <- noisy_boot_means_scaled[,i]
  q1 <- quantile(curr_data, 0.25)
  q2 <- quantile(curr_data, 0.5)
  q3 <- quantile(curr_data, 0.75)
  grid_seq <- seq(from=q1-3*(q3-q1), to=q3+3*(q3-q1), length.out = bin_num)
  
  clean_density <- density(clean_boot_means[,i])
  noisy_density <- density(curr_data * sd_of_noise)
  
  ### deconvolveR
  if (test_method == 'deconvolveR'){
    deconv_result <- deconv(tau = grid_seq, X = curr_data, family = "Normal", pDegree = 5, c0=c0)
    x_seq <- deconv_result$stats[, 'theta']
    CDF_seq <- deconv_result$stats[, 'G']
    PDF_seq <- deconv_result$stats[, 'g']
    max(PDF_seq) / PDF_seq[1]
  }
  
  den_x <- c(den_x, clean_density$x, noisy_density$x, x_seq * sd_of_noise)
  den_y <- c(den_y, clean_density$y, noisy_density$y, PDF_seq / sum(PDF_seq) / (
    (x_seq[3] - x_seq[2]) * sd_of_noise))
  facet_idx <- c(facet_idx, rep(i, length(clean_density$x) + length(noisy_density$x) + length(x_seq)))
  empirical_density_color <- c(empirical_density_color, c(rep("non-private bootstrap", length(clean_density$x)), 
                            rep("private bootstrap", length(clean_density$x)), 
                            rep("deconvolved private bootstrap", length(x_seq))))
  empirical_density_linetype <- c(empirical_density_linetype, c(rep(1, length(clean_density$x)), 
                                                                rep(2, length(noisy_density$x)), 
                                                                rep(3, length(x_seq))))
}
empirical_density_color <- factor(empirical_density_color, levels = c("non-private bootstrap",
                                                                      "private bootstrap",
                                                                      "deconvolved private bootstrap"))
empirical_density_linetype <- factor(empirical_density_linetype, levels=c(1,2,3))
deconv_df <- data.frame(den_x, den_y, empirical_density_color, 
                        empirical_density_linetype, facet_idx)

alpha_val <- 0.9
override.alpha <- c(alpha_val, alpha_val, alpha_val)
override.linetype <- c(1,3,2)
ggplot(data=deconv_df, aes(group=empirical_density_color, color=empirical_density_color, x=den_x, y=den_y)) + 
    geom_line(alpha=alpha_val, aes(linetype=empirical_density_linetype)) + xlim(0.48,0.52) + ylim(0,160) +
    theme_classic() + theme(legend.position="right", panel.background = element_rect(fill='transparent'),
                            legend.background = element_rect(fill='transparent'),
                            plot.background = element_rect(fill='transparent', color=NA),
                            legend.box.margin = margin(t=-0.0, unit='cm')) +
    xlab('theta') + ylab('density')  + 
    scale_linetype_manual(values=c('solid','dotted','dashed','dotdash','longdash'))+ # ,'twodash'
    facet_wrap(~ facet_idx, nrow = 1) +  
    scale_x_continuous(breaks = c(0.48,0.5,0.52)) +
    guides(color=guide_legend(title="Sampling distribution estimates", ncol=1, 
                              override.aes = list(alpha = override.alpha, linetype = override.linetype)), 
          linetype="none", alpha="none") 
ggsave(paste("Figure4.pdf", sep=""), width=7.7, height=1.5)







