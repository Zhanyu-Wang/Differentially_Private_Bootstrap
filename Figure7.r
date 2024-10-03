library(ggplot2)
library(ggpubr)
library(reshape2)
library(cowplot) 
library(tidyverse)



library(tidyverse)

nSIM <- 2000
n <- 10000  # sample size
# B <- 100 # 
ep_list <- c(1, 0.5, 0.3, 0.1) 

linetype_values <- c('solid','dashed','dotted','dotdash','longdash')
override.linetype <- c(1:length(linetype_values))

width_list = c()
for(ep in ep_list){
  B = as.integer(ep * ep * 2000)
  filename2 = paste("results/simulation_mean",
                    "_N=", n, "_mu=", ep, "_B=", B, "_nSIM=", nSIM, ".txt", sep='')
  if (!file.exists(filename2)) {
    print(filename2)
    stop()
  }
  CIs = read.csv(filename2)
  if (ep == ep_list[1]){
    width_list = rbind(width_list, c(as.double(CIs[1,2:5]), ep, "Non-private Bootstrap"))
    bootstrap_result = as.double(CIs[1,2:5])
  }
  else{
    width_list = rbind(width_list, c(bootstrap_result, ep, "Non-private Bootstrap"))
  }
  width_list = rbind(width_list, c(as.double(CIs[2,2:5]), ep, paste("DP Bootstrap (deconvolution)", sep="")))
  width_list = rbind(width_list, c(as.double(CIs[3,2:5]), ep, paste("DP Bootstrap (asymptotics)", sep="")))
  width_list = rbind(width_list, c(as.double(CIs[4,2:5]), ep, paste("Brawner and Honaker (2018)", sep="")))
  width_list = rbind(width_list, c(as.double(CIs[5,2:5]), ep, paste("NoisyVar (Du et al., 2020)", sep="")))
}
width_list = as.data.frame(width_list)
colnames(width_list) = c("coverage_0.9", "coverage_std_0.9", "mean_width_0.9", "std_width_0.9",
                         "ep","method")
width_list[, 1:(ncol(width_list)-1)] <- sapply(width_list[, 1:(ncol(width_list)-1)], as.double)
# width_list$ep = as.integer(width_list$ep)
width_list$method = factor(as.character(width_list$method), level=c("Non-private Bootstrap", 
                                                                    "DP Bootstrap (deconvolution)", 
                                                                    "DP Bootstrap (asymptotics)", 
                                                                    "Brawner and Honaker (2018)", 
                                                                    "NoisyVar (Du et al., 2020)"))
width_list$err_width_0.9 = 2*width_list$std_width_0.9
width_list$err_coverage_0.9 = 2*sqrt((1-width_list$coverage_0.9)*width_list$coverage_0.9/nSIM)


##############################
p1 <- ggplot(data=width_list, aes(x = ep, y = mean_width_0.9, group=method, lt=method, colour=(method))) +
  geom_line(aes(linetype=method)) +
  geom_errorbar(aes(ymin=pmax(0.0002, mean_width_0.9-err_width_0.9), ymax=mean_width_0.9+err_width_0.9), width=.01,
                position=position_dodge(.00)) +
  scale_y_continuous(limits = c(0.01, 0.033), breaks = c(0.01,0.02,0.03)) +
  scale_x_continuous(limits = c(0.08, 1.05), breaks = sort(ep_list)) +
  scale_color_manual(values=c('darkgoldenrod1','darkcyan','darkgreen','red','purple')) +
  labs(title=" CI width ",
       x ="GDP", y = "Population mean") +
  scale_linetype_manual(values=c('solid','dashed','dotted','dotdash','longdash'))+ # ,'twodash'
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), legend.position="top") + 
  guides(color=guide_legend(reverse = FALSE, title="Method", ncol=1, 
                            override.aes = list(linetype = override.linetype)), 
         linetype="none", alpha="none") 
p1

p2 <- ggplot(data=width_list, aes(x = ep, y = coverage_0.9, group=method, lt=method, colour=(method))) +
  geom_line(aes(linetype=method)) +
  geom_errorbar(aes(ymin=pmax(0, coverage_0.9-err_coverage_0.9), ymax=coverage_0.9+err_coverage_0.9), width=.01,
                position=position_dodge(.03)) +
  scale_y_continuous(limits = c(0.78, 1), breaks = c(0.8,0.9,1)) +
  scale_x_continuous(limits = c(0.08, 1.05), breaks = sort(ep_list)) +
  scale_color_manual(values=c('darkgoldenrod1','darkcyan','darkgreen','red','purple')) +
  labs(title="CI coverage",
       x ="GDP", y = NULL) +
  scale_linetype_manual(values=c('solid','dashed','dotted','dotdash','longdash'))+ # ,'twodash'
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), legend.position="top", ) + 
  guides(color=guide_legend(reverse = FALSE, title=expression(epsilon*"  "), nrow=1, 
                            override.aes = list(linetype = override.linetype)), 
         linetype="none", alpha="none") 



sp1 <- ggarrange(p1, p2, widths = c(3.8,3.4), ncol=2, nrow=1, common.legend = TRUE, legend="right")
sp1

figure_width <- 7
figure_height <- 2.2
ggsave("Figure7.pdf", width=figure_width, height=figure_height)










