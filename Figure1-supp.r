library(ggplot2)
library(tidyverse)
####################################################################

nSIM <- 2000
penalty_c <- 1
n_list <- c(100,300,1000,3000,10000,30000,100000)  # sample size
B_list <- c(1000,500,200,100,50,30,20,10,5) # 
ep_list <- c(0.1,0.2,0.5,1) #,1e9
save_str <- "clamp_normal_mean_comparison"
figure_width <- 14
figure_height <- 5
conf_level <- 0.9

linetype_values <- c('solid','dashed','dotted','dotdash','longdash','twodash','solid','dashed','dotted','dotdash','longdash')
override.linetype <- c(1:length(linetype_values))

width_list = c()
for(n in n_list){
  for(ep in ep_list){
    for(B in B_list){
      filename2 = paste("results_clamp_normal_mean/deconvolution_clamp_normal_mean_N(0.5,1)_[0,1]_", 
                        "_N=", n, "_mu=", ep, "_B=", B, "_nSIM=", nSIM, "_", conf_level, ".csv", sep='')
      if (!file.exists(filename2)) {
        print(filename2)
        stop()
      }
      CIs = read.csv(filename2)
      add_base_noise = 1 / n / ep
      add_dpboot_noise = 1 / n * 1.125 / ep * sqrt(B)
      if (B == B_list[1]){
        clean_width = CIs[1,4]
        clean_boot_std = clean_width / 2 / qnorm((1+conf_level)/2)
        width_list = rbind(width_list, c(as.double(CIs[3,2:5]), n, ep, (CIs[3,4]-clean_width)/clean_width, 1, "NoisyVar"))
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
width_list$B = factor(as.character(width_list$B), levels = c("Bootstrap (B=1000)","B=5","B=10","B=20","B=30","B=50","B=100","B=200","B=500","B=1000","NoisyVar"))
width_list$err_width = 2*width_list$std_width
width_list$err_coverage = 2*sqrt((1-width_list$coverage)*width_list$coverage/nSIM)
width_list$err_rej_rate = 2*sqrt((1-width_list$rej_rate)*width_list$rej_rate/nSIM)


privacy_labs <- c("mu=0.1","mu=0.2","mu=0.5","mu=1")
names(privacy_labs) <- c("0.1", "0.2", "0.5", "1")

##############################
p1 <- ggplot(data=width_list, aes(x = n, y = mean_width, group=B, lt=B, colour=factor(B))) +
  geom_line(aes(linetype=B)) +
  geom_errorbar(aes(ymin=pmax(0.003, mean_width-err_width), ymax=mean_width+err_width), width=.01,
                position=position_dodge(.01)) +
  scale_y_log10(limits = c(0.003, 6), breaks = c(0.01,0.1,1), labels = function(x) format(x, scientific = TRUE)) +
  scale_x_log10(limits = c(90, 120000), breaks = n_list) +
  labs(title=NULL,
       x ="sample size", y = paste("90% CI Width")) +
  scale_linetype_manual(values=linetype_values)+
  theme_bw() + theme(plot.title = element_text(hjust = 0), legend.position="top",axis.title.x=element_blank()) + 
  guides(color=guide_legend(reverse = FALSE, title=expression("  "), nrow=1, 
                            override.aes = list(linetype = override.linetype)), 
         linetype="none", alpha="none") +  
  facet_wrap(.~mu, nrow=1, labeller=labeller(mu = privacy_labs))

p2 <- ggplot(data=width_list, aes(x = n, y = coverage, group=B, lt=B, colour=factor(B))) +
  geom_line(aes(linetype=B)) +
  geom_errorbar(aes(ymin=pmax(0, coverage-err_coverage), ymax=coverage+err_coverage), width=.01,
                position=position_dodge(.07)) +
  scale_y_continuous(breaks = c(0.8,0.9,1), labels = function(x) format(x, scientific = TRUE)) +
  scale_x_log10(limits = c(90, 120000), breaks = n_list) +
  labs(title=NULL,
       x ="sample size", y = paste("90% CI Coverage")) +
  scale_linetype_manual(values=linetype_values)+
  theme_bw() + theme(plot.title = element_text(hjust = 0), legend.position="top",axis.title.x=element_blank()) + 
  guides(color=guide_legend(reverse = FALSE, title=expression("  "), nrow=1, 
                            override.aes = list(linetype = override.linetype)), 
         linetype="none", alpha="none") +  
  facet_wrap(.~mu, nrow=1, labeller=labeller(mu = privacy_labs))

p3 <- ggplot(data=width_list, aes(x = n, y = pmax(0.00001,width_ratio), group=B, lt=B, colour=factor(B))) +
  geom_line(aes(linetype=B)) +
  scale_y_log10(limits = c(0.00001, 30), breaks = c(0.0001,0.001,0.01,0.1,1,10), labels = function(x) format(x, scientific = TRUE)) +
  scale_x_log10(limits = c(90, 120000), breaks = n_list) +
  labs(title=NULL,
       x ="sample size", y = "Extra width") +
  scale_linetype_manual(values=linetype_values)+
  theme_bw() + theme(plot.title = element_text(hjust = 0), legend.position="top", axis.title.x=element_blank()) + 
  guides(color=guide_legend(reverse = FALSE, title=expression("  "), nrow=1, 
                            override.aes = list(linetype = override.linetype)), 
         linetype="none", alpha="none") +  
  facet_wrap(.~mu, nrow=1, labeller=labeller(mu = privacy_labs))

p4 <- ggplot(data=width_list, aes(x = n, y = noise_ratio, group=B, lt=B, colour=factor(B))) +
  geom_line(aes(linetype=B)) +
  geom_hline(yintercept=1) +
  scale_y_log10(limits = c(0.01, 100), breaks = c(0.01,0.1,1,10,100), labels = function(x) format(x, scientific = TRUE)) +
  scale_x_log10(limits = c(90, 120000), breaks = n_list) +
  labs(title=NULL,
       x ="sample size", y = "sqrt(SNR)") +
  scale_linetype_manual(values=linetype_values)+
  theme_bw() + theme(plot.title = element_text(hjust = 0), legend.position="top") + 
  guides(color=guide_legend(reverse = FALSE, title=expression("  "), nrow=1, 
                            override.aes = list(linetype = override.linetype)), 
         linetype="none", alpha="none") +  
  facet_wrap(.~mu, nrow=1, labeller=labeller(mu = privacy_labs))



library(ggpubr)
ggarrange(p1, p2, p3, p4, ncol=1, nrow=4, common.legend = TRUE, legend="top")
figure_width <- 12
figure_height <- 7
ggsave(paste(save_str, "_", conf_level, ".pdf", sep=""), width=figure_width, height=figure_height)








