library(ggplot2)
library(ggpubr)
library(reshape2)
library(cowplot) 
library(tidyverse)

### 'ind2016.csv' is the original file we download from https://www.kaggle.com/datasets/smitop/statcan-pumf-2016

# ind_data <- read.csv(file = 'ind2016.csv') 
# ind_data2 <- ind_data[,c(94,115,130)]
# write.csv(ind_data2, 'ind2016_MRKINC-PR-SHELCO.csv')
ind_data2 <- read.csv(file = 'ind2016_MRKINC-PR-SHELCO.csv')

description_str <- "ind2016-cov_"

row_col <- dim(ind_data2) # 930421    4

# col_index:  MRKINC=94, SHELCO=130
# PR:         Alberta=48, Quebec=24, Ontario=35, British Columbia=59.
# provide_PR_codes <- data.frame(Alberta=48, Quebec=24, Ontario=35, British_Columbia=59)
provide_PR_codes <- data.frame(Quebec=24, Ontario=35)
col_indices <- data.frame(MRKINC=94, SHELCO=130)
# col_indices <- data.frame(MRKINC=94)

provide_PR_code_name = 'Ontario'

provide_PR_code <- as.numeric(provide_PR_codes[provide_PR_code_name])
province_data <- ind_data2[ind_data2[1:row_col[1], "PR"] == provide_PR_code, colnames(col_indices)]
row_indice <- (province_data != 99999999) & (province_data != 88888888) &
  (province_data != 678800) & (province_data %% 100 == 0) & (province_data >= 0)
row_indice_1 <- row_indice[,1] & row_indice[,2]
province_data <- province_data[row_indice_1,]
hull <- province_data %>%
  slice(chull(MRKINC, SHELCO))
sp2 <- ggplot(province_data, aes(x = MRKINC, y = SHELCO)) +
  geom_polygon(data = hull, alpha = 1, fill=('lightblue')) +
  stat_density_2d(aes(fill = (..level..)), bins = 100, geom = "polygon") +
  gradient_fill(c("lightblue", "steelblue")) + guides(fill="none") +
  theme_bw() +
  scale_x_continuous(breaks = c(0, max(province_data$MRKINC)), labels = scales::comma) +
  scale_y_continuous(breaks = c(0, max(province_data$SHELCO)), labels = scales::comma) +
  theme(legend.position="top", panel.background = element_rect(fill='transparent'),
        legend.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        plot.margin = margin(t = 10, r = 20, l = 10, b = 10, unit = "pt"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  coord_cartesian(xlim=c(0, max(province_data$MRKINC)), ylim=c(0, max(province_data$SHELCO)), expand=0) +
  ggtitle(NULL) # provide_PR_code_name





library(tidyverse)

nSIM <- 2000
penalty_c <- 1
n_list <- c(1000,3000,10000,30000,100000)  # sample size
B <- 100 # 
ep_list <- c(1) # 0.01, 0.03, 0.1, 0.3, 
description_str <- "ind2016-logistic_"
# provide_PR_code_name <- "Quebec"
provide_PR_code_name <- "Ontario"
figure_width <- 10
figure_height <- 4

linetype_values <- c('solid','dashed','dotted','dotdash')
override.linetype <- c(1:length(linetype_values))

width_list = c()
for(penalty_c in c(penalty_c)){
  for(n in n_list){
    for(ep in ep_list){
      filename2 = paste("results_logistic/deconvolution_", description_str, provide_PR_code_name, 
                        "_N=", n, "_mu=", ep, "_B=", B, "_c=", penalty_c, "_nSIM=", nSIM, ".csv", sep='')
      if (!file.exists(filename2)) {
        print(filename2)
        stop()
      }
      CIs = read.csv(filename2)
      if (ep == ep_list[1]){
        width_list = rbind(width_list, c(as.double(CIs[1,2:5]), n, "Bootstrap"))
      }
      width_list = rbind(width_list, c(as.double(CIs[4,2:5]), n, paste("DP-CI-ERM (Wang et al., 2019)", sep="")))
      width_list = rbind(width_list, c(as.double(CIs[3,2:5]), n, paste("DP Bootstrap (deconvolution)", sep="")))
      width_list = rbind(width_list, c(as.double(CIs[2,2:5]), n, paste("DP Bootstrap (asymptotics)", sep="")))
    }
  }
}
width_list = as.data.frame(width_list)
colnames(width_list) = c("coverage_0.9", "rej_rate_0.9", "mean_width_0.9", "std_width_0.9",
                         "n","method")
width_list[, 1:(ncol(width_list)-1)] <- sapply(width_list[, 1:(ncol(width_list)-1)], as.double)
width_list$n = as.integer(width_list$n)
width_list$method = factor(as.character(width_list$method), level=c("Bootstrap", "DP Bootstrap (asymptotics)", "DP Bootstrap (deconvolution)", "DP-CI-ERM (Wang et al., 2019)"))
width_list$err_width_0.9 = 2*width_list$std_width_0.9
width_list$err_coverage_0.9 = 2*sqrt((1-width_list$coverage_0.9)*width_list$coverage_0.9/nSIM)
width_list$err_rej_rate_0.9 = 2*sqrt((1-width_list$rej_rate_0.9)*width_list$rej_rate_0.9/nSIM)


##############################
p1 <- ggplot(data=width_list, aes(x = n, y = mean_width_0.9, group=method, lt=method, colour=(method))) +
  geom_line(aes(linetype=method)) +
  geom_errorbar(aes(ymin=pmax(0.0002, mean_width_0.9-err_width_0.9), ymax=mean_width_0.9+err_width_0.9), width=.01,
                position=position_dodge(.01)) +
  scale_y_log10(limits = c(0.0002, 0.03), breaks = c(0.0003,0.001,0.003,0.01,0.03)) +
  scale_x_log10(limits = c(900, 120000), breaks = n_list, labels = scales::label_number_si()) +
  scale_color_manual(values=c('darkgoldenrod1','darkcyan','darkgreen','red')) +
  labs(title=" CI width ",
       x ="sample size", y = "Logistic \n regression") +
  scale_linetype_manual(values=c('solid','dashed','dotted','dotdash','longdash'))+ # ,'twodash'
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), legend.position="top", 
        axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x=element_blank()) + 
  guides(color=guide_legend(reverse = FALSE, title="Method", ncol=1, 
                            override.aes = list(linetype = override.linetype)), 
         linetype="none", alpha="none") 


p2 <- ggplot(data=width_list, aes(x = n, y = coverage_0.9, group=method, lt=method, colour=(method))) +
  geom_line(aes(linetype=method)) +
  geom_errorbar(aes(ymin=pmax(0, coverage_0.9-err_coverage_0.9), ymax=coverage_0.9+err_coverage_0.9), width=.01,
                position=position_dodge(.07)) +
  scale_y_continuous(limits = c(0.87, 1), breaks = c(0.9,0.95,1)) +
  scale_x_log10(limits = c(900, 120000), breaks = n_list, labels = scales::label_number_si()) +
  scale_color_manual(values=c('darkgoldenrod1','darkcyan','darkgreen','red')) +
  labs(title="CI coverage",
       x ="sample size", y = NULL) +
  scale_linetype_manual(values=c('solid','dashed','dotted','dotdash','longdash'))+ # ,'twodash'
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), legend.position="top", 
        axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x=element_blank()) + 
  guides(color=guide_legend(reverse = FALSE, title=expression(epsilon*"  "), nrow=1, 
                            override.aes = list(linetype = override.linetype)), 
         linetype="none", alpha="none") 


p3 <- ggplot(data=width_list, aes(x = n, y = rej_rate_0.9, group=method, lt=method, colour=(method))) +
  geom_line(aes(linetype=method)) +
  geom_errorbar(aes(ymin=pmax(0, rej_rate_0.9-err_rej_rate_0.9), ymax=rej_rate_0.9+err_rej_rate_0.9), width=.01,
                position=position_dodge(.01)) +
  # scale_y_log10(limits = c(0.9, 1.1), breaks = c(0.002,0.01,0.05,0.2,1,5)) +
  scale_x_log10(limits = c(900, 120000), breaks = n_list, labels = scales::label_number_si()) +
  scale_color_manual(values=c('darkgoldenrod1','darkcyan','darkgreen','red')) +
  labs(title="P(CI covers 0)",
       x ="sample size", y = NULL) +
  scale_linetype_manual(values=c('solid','dashed','dotted','dotdash','longdash'))+ # ,'twodash'
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), legend.position="top", 
        axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x=element_blank()) + 
  guides(color=guide_legend(reverse = FALSE, title=expression(epsilon*"  "), nrow=1, 
                            override.aes = list(linetype = override.linetype)), 
         linetype="none", alpha="none") 


description_str <- "ind2016-quantile_"

linetype_values2 <- c('solid','dashed','dotted')
override.linetype2 <- c(1:length(linetype_values2))

width_list2 = c()
for(penalty_c in c(penalty_c)){
  for(B in c(B)){
    for(n in n_list){
      for(ep in ep_list){
        filename2 = paste("results_quantile/deconvolution_", description_str, provide_PR_code_name, 
                          "_N=", n, "_mu=", ep, "_B=", B, "_c=", penalty_c, "_nSIM=", nSIM, ".csv", sep='')
        if (!file.exists(filename2)) {
          print(filename2)
          stop()
        }
        CIs = read.csv(filename2)
        if (ep == ep_list[1]){
          width_list2 = rbind(width_list2, c(as.double(CIs[1,2:5]), n, penalty_c, B, "Bootstrap"))
        }
        width_list2 = rbind(width_list2, c(as.double(CIs[3,2:5]), n, penalty_c, B, paste("DP Bootstrap (deconvolution)",sep="")))
        width_list2 = rbind(width_list2, c(as.double(CIs[2,2:5]), n, penalty_c, B, paste("DP Bootstrap (asymptotics)",sep="")))
      }
    }
  }
}
width_list2 = as.data.frame(width_list2)
colnames(width_list2) = c("coverage_0.9", "rej_rate_0.9", "mean_width_0.9", "std_width_0.9",
                          "n","c","B","ep")
width_list2[, 1:(ncol(width_list2)-1)] <- sapply(width_list2[, 1:(ncol(width_list2)-1)], as.double)
width_list2$n = as.integer(width_list2$n)
width_list2$method = factor(as.character(width_list2$ep), levels = c("Bootstrap", "DP Bootstrap (asymptotics)", "DP Bootstrap (deconvolution)", "DP-CI-ERM"))
width_list2$err_width_0.9 = 2*width_list2$std_width_0.9
width_list2$err_coverage_0.9 = 2*sqrt((1-width_list2$coverage_0.9)*width_list2$coverage_0.9/nSIM)
width_list2$err_rej_rate_0.9 = 2*sqrt((1-width_list2$rej_rate_0.9)*width_list2$rej_rate_0.9/nSIM)


##############################
p4 <- ggplot(data=width_list2, aes(x = n, y = mean_width_0.9, group=method, lt=method, colour=(method))) +
  geom_line(aes(linetype=method)) +
  geom_errorbar(aes(ymin=pmax(0.0004, mean_width_0.9-err_width_0.9), ymax=mean_width_0.9+err_width_0.9), width=.01,
                position=position_dodge(.01)) +
  scale_y_log10(limits = c(0.0004, 0.0205), breaks = c(0.0005,0.002,0.01)) +
  scale_x_log10(limits = c(900, 120000), breaks = n_list, labels = scales::label_number_si()) +
  scale_color_manual(values=c('darkgoldenrod1','darkcyan','darkgreen','red')) +
  labs(title=NULL,
       x ="sample size", y = "Quantile \n regression") +
  scale_linetype_manual(values=c('solid','dashed','dotted','dotdash','longdash'))+ # ,'twodash'
  theme_bw() + theme(plot.title = element_text(hjust = 0), legend.position="top") + 
  guides(color=guide_legend(reverse = FALSE, title="method", nrow=1, 
                            override.aes = list(linetype = override.linetype2)), 
         linetype="none", alpha="none")


p5 <- ggplot(data=width_list2, aes(x = n, y = coverage_0.9, group=method, lt=method, colour=(method))) +
  geom_line(aes(linetype=method)) +
  geom_errorbar(aes(ymin=pmax(0, coverage_0.9-err_coverage_0.9), ymax=coverage_0.9+err_coverage_0.9), width=.01,
                position=position_dodge(.07)) +
  scale_y_continuous(limits = c(0.87, 1), breaks = c(0.9,0.95,1)) +
  scale_x_log10(limits = c(900, 120000), breaks = n_list, labels = scales::label_number_si()) +
  scale_color_manual(values=c('darkgoldenrod1','darkcyan','darkgreen','red')) +
  labs(title=NULL,
       x ="sample size", y = NULL) +
  scale_linetype_manual(values=c('solid','dashed','dotted','dotdash','longdash'))+ # ,'twodash'
  theme_bw() + theme(plot.title = element_text(hjust = 0), legend.position="top") + 
  guides(color=guide_legend(reverse = FALSE, title=expression(epsilon*"  "), nrow=1, 
                            override.aes = list(linetype = override.linetype2)), 
         linetype="none", alpha="none") 


p6 <- ggplot(data=width_list2, aes(x = n, y = rej_rate_0.9, group=method, lt=method, colour=(method))) +
  geom_line(aes(linetype=method)) +
  geom_errorbar(aes(ymin=pmax(0, rej_rate_0.9-err_rej_rate_0.9), ymax=rej_rate_0.9+err_rej_rate_0.9), width=.01,
                position=position_dodge(.01)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25,0.5,0.75,1)) +
  scale_x_log10(limits = c(900, 120000), breaks = n_list, labels = scales::label_number_si()) +
  scale_color_manual(values=c('darkgoldenrod1','darkcyan','darkgreen','red')) +
  labs(title=NULL,
       x ="sample size", y = NULL) +
  scale_linetype_manual(values=c('solid','dashed','dotted','dotdash','longdash'))+ # ,'twodash'
  theme_bw() + theme(plot.title = element_text(hjust = 0), legend.position="top") + 
  guides(color=guide_legend(reverse = FALSE, title=expression(epsilon*"  "), nrow=1, 
                            override.aes = list(linetype = override.linetype2)), 
         linetype="none", alpha="none") 


sp1 <- ggarrange(p1, p2, p3, p4, p5, p6, widths = c(3.8,3.,3.), ncol=3, nrow=2, common.legend = TRUE, legend="right")






ggarrange(sp2, sp1, ncol = 2, nrow = 1, widths = c(0.8,2.2), labels = "auto")

figure_width <- 11
figure_height <- 2.8
ggsave("Figure8.pdf", width=figure_width, height=figure_height)










