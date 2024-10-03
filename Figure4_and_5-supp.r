library(ggplot2)
library(tidyverse)
####################################################################

nSIM <- 2000
penalty_c <- 1
n_list <- c(1000,3000,10000,30000,100000)  # sample size
B <- 200 # 
ep_list <- c(1,2,5,10) 
B_list <- c(1000,500,200,100,50,20)
description_str <- "ind2016-quantile_"
save_str <- "quantile_regression_comparison"
provide_PR_code_name <- "Ontario"
figure_width <- 10
figure_height <- 10
conf_level <- 0.9

linetype_values <- c('solid','dashed','dotted','dotdash','longdash','twodash','solid')
override.linetype <- c(1:length(linetype_values))
cbbPalette <- c("#F8766C", "#CD9500", "#7BAD00", "#00BD67", "#00BEC4", "#00A9FF", "#C67BFF", "#FF60CC")

for(penalty_c in c(1,0.01)){
  width_list = c()
  for(n in n_list){
    for(ep in ep_list){
      for(B in B_list){
        filename2 = paste("results_quantile/deconvolution_", description_str, provide_PR_code_name, 
                          "_N=", n, "_mu=", ep, "_B=", B, "_c=", penalty_c, "_nSIM=", nSIM, ".csv", sep='')
        if (!file.exists(filename2)) {
          print(filename2)
          stop()
        }
        CIs = read.csv(filename2)
        add_base_noise = 1 / n / (ep/2) / penalty_c
        add_dpboot_noise = 1 / n * 1.125 / ep * sqrt(B) / penalty_c
        if (B == B_list[1]){
          clean_width = CIs[1,4]
          clean_boot_std = clean_width / 2 / qnorm((1+conf_level)/2)
          width_list = rbind(width_list, c(as.double(CIs[1,2:5]), n, ep, (CIs[1,4]-clean_width)/clean_width, 0, "Bootstrap (B=1000)"))
        }
        width_list = rbind(width_list, c(as.double(CIs[2,2:5]), n, ep, (CIs[2,4]-clean_width)/clean_width, clean_boot_std / add_dpboot_noise, paste("B=",B,sep="")))
      }
    }
  }
  width_list = as.data.frame(width_list)
  colnames(width_list) = c("coverage", "rej_rate", "mean_width", "std_width",
                           "n","mu","width_ratio","noise_ratio","B")
  width_list[, 1:(ncol(width_list)-1)] <- sapply(width_list[, 1:(ncol(width_list)-1)], as.double)
  width_list$n = as.integer(width_list$n)
  width_list$B = factor(as.character(width_list$B), levels = c("Bootstrap (B=1000)","B=20","B=50","B=100","B=200","B=500","B=1000","DP-CI-ERM"))
  width_list$err_width = 2*width_list$std_width
  width_list$err_coverage = 2*sqrt((1-width_list$coverage)*width_list$coverage/nSIM)
  width_list$err_rej_rate = 2*sqrt((1-width_list$rej_rate)*width_list$rej_rate/nSIM)
  

  privacy_labs <- c("mu=1","mu=2","mu=5","mu=10")
  names(privacy_labs) <- c("1", "2", "5", "10")

  ##############################
  p1 <- ggplot(data=width_list, aes(x = n, y = mean_width, group=B, lt=B, colour=(B))) +
    geom_line(aes(linetype=B)) +
    geom_errorbar(aes(ymin=pmax(0.0001, mean_width-err_width), ymax=mean_width+err_width), width=.01,
                  position=position_dodge(.01)) +
    scale_y_log10(limits = c(0.0001, 10), breaks = c(0.0001,0.001,0.01,0.1,1), labels = function(x) format(x, scientific = TRUE)) +
    scale_x_log10(limits = c(800, 120000), breaks = n_list) +
    scale_colour_manual(values=cbbPalette) +
    labs(title=NULL,
         x ="sample size", y = paste("90% CI width")) +
    scale_linetype_manual(values=linetype_values)+
    theme_bw() + theme(plot.title = element_text(hjust = 0), legend.position="top",axis.title.x=element_blank()) + 
    guides(color=guide_legend(reverse = FALSE, title=expression("  "), nrow=1, 
                              override.aes = list(linetype = override.linetype)), 
           linetype="none", alpha="none") +  
    facet_wrap(.~mu, nrow=1, labeller=labeller(mu = privacy_labs))
  
  p2 <- ggplot(data=width_list, aes(x = n, y = coverage, group=B, lt=B, colour=(B))) +
    geom_line(aes(linetype=B)) +
    geom_errorbar(aes(ymin=pmax(0, coverage-err_coverage), ymax=coverage+err_coverage), width=.01,
                  position=position_dodge(.07)) +
    scale_y_continuous(breaks = c(0.8,0.9,1), labels = function(x) format(x, scientific = TRUE)) +
    scale_x_log10(limits = c(800, 120000), breaks = n_list) +
    scale_colour_manual(values=cbbPalette) +
    labs(title=NULL,
         x ="sample size", y = paste("90% CI Coverage")) +
    scale_linetype_manual(values=linetype_values)+
    theme_bw() + theme(plot.title = element_text(hjust = 0), legend.position="top",axis.title.x=element_blank()) + 
    guides(color=guide_legend(reverse = FALSE, title=expression("  "), nrow=1, 
                              override.aes = list(linetype = override.linetype)), 
           linetype="none", alpha="none") +  
    facet_wrap(.~mu, nrow=1, labeller=labeller(mu = privacy_labs))
    
  p3 <- ggplot(data=width_list, aes(x = n, y = rej_rate, group=B, lt=B, colour=(B))) +
    geom_line(aes(linetype=B)) +
    geom_errorbar(aes(ymin=pmax(0, rej_rate-err_rej_rate), ymax=rej_rate+err_rej_rate), width=.01,
                  position=position_dodge(.01)) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0.,0.5,1), labels = function(x) format(x, scientific = TRUE)) +
    scale_x_log10(limits = c(800, 120000), breaks = n_list) +
    scale_colour_manual(values=cbbPalette) +
    labs(title=NULL,
         x ="sample size", y = "P(CI covers 0)") +
    scale_linetype_manual(values=linetype_values) +
    theme_bw() + theme(plot.title = element_text(hjust = 0), legend.position="top",axis.title.x=element_blank()) + 
    guides(color=guide_legend(reverse = FALSE, title=expression("  "), nrow=1, 
                              override.aes = list(linetype = override.linetype)), 
           linetype="none", alpha="none") +  
    facet_wrap(.~mu, nrow=1, labeller=labeller(mu = privacy_labs))
  
  p4 <- ggplot(data=width_list, aes(x = n, y = pmax(1e-4,width_ratio), group=B, lt=B, colour=(B))) +
    geom_line(aes(linetype=B)) +
    scale_y_log10(limits = c(1e-4, 30), breaks = c(0.0001,0.001,0.01,0.1,1,10), labels = function(x) format(x, scientific = TRUE)) +
    scale_x_log10(limits = c(800, 120000), breaks = n_list) +
    scale_colour_manual(values=cbbPalette) +
    labs(title=NULL,
         x ="sample size", y = paste("Extra width")) +
    scale_linetype_manual(values=linetype_values)+
    theme_bw() + theme(plot.title = element_text(hjust = 0), legend.position="top",axis.title.x=element_blank()) + 
    guides(color=guide_legend(reverse = FALSE, title=expression("  "), nrow=1, 
                              override.aes = list(linetype = override.linetype)), 
           linetype="none", alpha="none") +  
    facet_wrap(.~mu, nrow=1, labeller=labeller(mu = privacy_labs))
  
  p5 <- ggplot(data=width_list, aes(x = n, y = noise_ratio, group=B, lt=B, colour=(B))) +
    geom_line(aes(linetype=B)) +
    geom_hline(yintercept=1) +
    scale_y_log10(limits = c(0.01, 30), breaks = c(0.01,0.01,0.1,1,10), labels = function(x) format(x, scientific = TRUE)) +
    scale_x_log10(limits = c(800, 120000), breaks = n_list) +
    scale_colour_manual(values=cbbPalette) +
    labs(title=NULL,
         x ="sample size", y = paste("sqrt(SNR)")) +
    scale_linetype_manual(values=linetype_values)+
    theme_bw() + theme(plot.title = element_text(hjust = 0), legend.position="top") + 
    guides(color=guide_legend(reverse = FALSE, title=expression("  "), nrow=1, 
                              override.aes = list(linetype = override.linetype)), 
           linetype="none", alpha="none") +  
    facet_wrap(.~mu, nrow=1, labeller=labeller(mu = privacy_labs))
  
  
  
  library(ggpubr)
  ggarrange(p1, p2, p3, p4, p5, ncol=1, nrow=5, common.legend = TRUE, legend="top")
  figure_width <- 12
  figure_height <- 8
  ggsave(paste(save_str, "_", conf_level, "_", provide_PR_code_name, "_", penalty_c, ".pdf", sep=""), width=figure_width, height=figure_height)

}











