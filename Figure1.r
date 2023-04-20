library(ggplot2)


total_len <- 1001
x_seq <- seq(-3,3,length=total_len)
typeI_seq <- pnorm(-x_seq, mean=0.5)
typeII_seq <- pnorm(x_seq, mean=0.5)
typeI_seq2 <- pnorm(-x_seq, mean=0.65)
typeII_seq2 <- pnorm(x_seq, mean=0.65)
typeI_seq3_1 <- 0.2+0.2/3*seq(-3,0,length=(total_len-1)/2)
typeII_seq3_1 <- 0.2-0.8/3*seq(-3,0,length=(total_len-1)/2)
typeI_seq3 <- c(typeI_seq3_1,0.2,rev(typeII_seq3_1))
typeII_seq3 <- c(typeII_seq3_1,0.2,rev(typeI_seq3_1))
typeI_seq4 <- pnorm(-x_seq, mean=0.3)
typeII_seq4 <- pnorm(x_seq, mean=0.3)
tradeoff.df <- data.frame(typeI_seq=c(typeI_seq, typeI_seq2,typeI_seq3,typeI_seq4), 
                          typeII_seq=c(typeII_seq, typeII_seq2,typeII_seq3,typeII_seq4),
                          privacy=rep(c("f", "Not f-DP (1)", "Not f-DP (2)", "Is f-DP"), each=length(typeI_seq)),
                          color=rep(c("red", "blue", "green", "yellow"), each=length(typeI_seq)))

map_str_to_int <- function(deconvolution_method){
  deconvolution_method[deconvolution_method == "f"] = 'solid'
  deconvolution_method[deconvolution_method == "Not f-DP (1)"] = 'dotted'
  deconvolution_method[deconvolution_method == "Not f-DP (2)"] = 'dashed'
  deconvolution_method[deconvolution_method == "Is f-DP"] = 'dotdashed'
  return(deconvolution_method)
  # return(as.character((deconvolution_method == "Not f-DP (1)") * 1 + (deconvolution_method == "Not f-DP (2)") * 10 + 
                        # (deconvolution_method == "Is f-DP") * 3))
}
override.alpha <- c(1,1,1,1)
override.linetype <- c(1,4,2,3)
ggplot() + 
  geom_line(data=tradeoff.df, mapping=aes(x=typeI_seq, y=typeII_seq, color=privacy, linetype=map_str_to_int(privacy))) +
  scale_linetype_manual(values = c("solid"=1, "dotted"=2, "dashed"=3, "dotdashed"=4)) +
  # scale_color_manual(values = c("red"="red", "blue"="blue", "green"="green", "purple"="purple")) +
  theme_classic() + theme(legend.position="right", panel.background = element_rect(fill='transparent'),
                          legend.background = element_rect(fill='transparent'),
                          plot.background = element_rect(fill='transparent', color=NA)) + 
  xlab('Type I Error') + ylab('Type II Error') +
  guides(color=guide_legend(title="Tradeoff functions", ncol=1,
                            override.aes = list(alpha = override.alpha, linetype = override.linetype)),
         linetype="none",
         alpha="none")
# scale_colour_manual(values = c("#000000", "#AC8254", "#5AA800", "#0098E9"), labels = c("f", "Is f-DP", "Not f-DP", "Not f-DP"))
ggsave("Figure1.pdf", width=3.6, height=2.0)
