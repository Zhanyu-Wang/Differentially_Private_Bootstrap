library(ggplot2)
library(latex2exp)
width = 3
i=0
nSIM <- 2000  # simulation number
N <- 10000  # sample size
B_0 <- 1000
pDegree <- 5
c0 <- 1e-1
bin_num <- 1000
readfile_all = c()
for (gdp_mu in c(1, 0.5, 0.3, 0.1)){
  i=i+1
  for (B in c(20, 180, 500, 2000)) {
    method_name <- paste('deconvolveR', sep='')
    csv_name <- paste("results/deconvolution_simulation-", method_name, "_range=1", 
                      "_N=", N, "_mu=", gdp_mu, "_B=", B,
                      "_pd=", pDegree, "_c0=", c0, "_bin=", bin_num, 
                      "_nSIM=", nSIM, ".csv", sep='')
    readfile <- read.csv(csv_name)
    readfile['facet_idx.labs'] <- rep(paste("mu=",gdp_mu,sep=""), length(readfile$variable))
    readfile['facet_idx'] <- rep(i, length(readfile))
    if(B == 2000){
      list_temp <- rep(paste("deconvolveR",B,sep='_'), length(readfile$variable))
      list_temp[readfile['deconvolution_method'] == 'clean_bootstrap'] <- rep("clean_bootstrap", length(readfile$variable)/2)
      readfile['method'] <- list_temp
      readfile_all <- rbind(readfile_all, readfile)
    } else {
      readfile['method'] <- rep(paste("deconvolveR",B,sep='_'), length(readfile))
      readfile_all <- rbind(readfile_all, readfile[readfile['deconvolution_method'] != 'clean_bootstrap',])
    }
  }
}

facet_idx.labs <- c("mu=1", "mu=0.5", "mu=0.3", "mu=0.1")
names(facet_idx.labs) <- c(1:length(facet_idx.labs)) 

linetype_values <- c('solid','dashed','dotted','dotdash','longdash')
override.linetype <- c(1:length(linetype_values))

readfile_all$method <- factor(as.character(readfile_all$method), levels = c("clean_bootstrap", "deconvolveR_2000", "deconvolveR_500", "deconvolveR_180", "deconvolveR_20"))
ggplot() + 
  stat_ecdf(alpha=1, data=readfile_all, aes(linetype=method, x=value, colour=method), geom = "step") + 
  theme_classic() + theme(legend.position="right", panel.background = element_rect(fill='transparent'),
                          legend.background = element_rect(fill='transparent'),
                          plot.background = element_rect(fill='transparent', color=NA)) +
  scale_colour_discrete(
    labels = c("clean_bootstrap" = "non-private bootstrap (B=2000)",
               "deconvolveR_20" = "deconvolved private bootstrap (B=20)",
               "deconvolveR_180" = "deconvolved private bootstrap (B=180)",
               "deconvolveR_500" = "deconvolved private bootstrap (B=500)",
               "deconvolveR_2000" = "deconvolved private bootstrap (B=2000)")
  ) +
  scale_linetype_manual(values=linetype_values) +
  facet_wrap(~ facet_idx, nrow = 1, labeller = labeller(facet_idx = facet_idx.labs)) +
  xlab("x=F*(theta)") + 
  ylab('F(x)')  +  
  scale_x_continuous(breaks = c(0.0,0.5,1.0)) +
  guides(color=guide_legend(title="Inference methods", ncol=1, 
                            override.aes = list(linetype = override.linetype)), 
         linetype="none", alpha="none") 
ggsave(paste("Figure5.pdf", sep=""), width=9, height=1.8)

