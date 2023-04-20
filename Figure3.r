library(rmutil)
library(ggplot2)
mu=1

N=100
k = seq(0,N)
p = dbinom(k,size=N,prob=1/N)

T = seq(-50,50, length=10001)

 
plot_width = 3
plot_height = 2.5
alpha_val = 0.6
line_type_val = 1

# 0 vs 1
shift = 0
typeI = function(t){
  k=seq(0,N)
  error=p[1:(N+1)]%*%(1-pnorm(t+k*mu*(shift)))
  return(error)
}
typeII = function(t){
  k=seq(0,N)
  error = p[1:(N+1)]%*%pnorm(t-k*mu*(1-shift))
  return(error)
}
Ivec1 = sapply(T,typeI)
IIvec1 = sapply(T,typeII)

mix_delta.df <- data.frame(Ivec1, IIvec1, line_type=1, alpha_val=alpha_val, data_pair=rep('0', length(IIvec1)))


# -0.2
shift = 0.2
typeI = function(t){
  k=seq(0,N)
  error=p[1:(N+1)]%*%(1-pnorm(t+k*mu*(shift)))
  return(error)
}
typeII = function(t){
  k=seq(0,N)
  error = p[1:(N+1)]%*%pnorm(t-k*mu*(1-shift))
  return(error)
}
Ivec2 = sapply(T,typeI)
IIvec2 = sapply(T,typeII)

mix_delta_new.df <- data.frame(Ivec1=Ivec2, IIvec1=IIvec2, line_type=line_type_val, alpha_val=alpha_val, data_pair=rep('1', length(IIvec1)))
mix_delta.df <- rbind(mix_delta.df, mix_delta_new.df)

# -0.4
shift = 0.4
typeI = function(t){
  k=seq(0,N)
  error=p[1:(N+1)]%*%(1-pnorm(t+k*mu*(shift)))
  return(error)
}
typeII = function(t){
  k=seq(0,N)
  error = p[1:(N+1)]%*%pnorm(t-k*mu*(1-shift))
  return(error)
}
Ivec2 = sapply(T,typeI)
IIvec2 = sapply(T,typeII)
mix_delta_new.df <- data.frame(Ivec1=Ivec2, IIvec1=IIvec2, line_type=line_type_val, alpha_val=rep(alpha_val, length(IIvec1)), data_pair=rep('2', length(IIvec1)))
mix_delta.df <- rbind(mix_delta.df, mix_delta_new.df)

# -0.4
shift = 0.5
typeI = function(t){
  k=seq(0,N)
  error=p[1:(N+1)]%*%(1-pnorm(t+k*mu*(shift)))
  return(error)
}
typeII = function(t){
  k=seq(0,N)
  error = p[1:(N+1)]%*%pnorm(t-k*mu*(1-shift))
  return(error)
}
Ivec2 = sapply(T,typeI)
IIvec2 = sapply(T,typeII)
mix_delta_new.df <- data.frame(Ivec1=Ivec2, IIvec1=IIvec2, line_type=4, alpha_val=rep(alpha_val, length(IIvec1)), data_pair=rep('3', length(IIvec1)))
mix_delta.df <- rbind(mix_delta.df, mix_delta_new.df)

# -0.6
shift = 0.6
typeI = function(t){
  k=seq(0,N)
  error=p[1:(N+1)]%*%(1-pnorm(t+k*mu*(shift)))
  return(error)
}
typeII = function(t){
  k=seq(0,N)
  error = p[1:(N+1)]%*%pnorm(t-k*mu*(1-shift))
  return(error)
}
Ivec2 = sapply(T,typeI)
IIvec2 = sapply(T,typeII)
mix_delta_new.df <- data.frame(Ivec1=Ivec2, IIvec1=IIvec2, line_type=line_type_val, alpha_val=rep(alpha_val, length(IIvec1)), data_pair=rep('4', length(IIvec1)))
mix_delta.df <- rbind(mix_delta.df, mix_delta_new.df)


# -0.8
shift = 0.8
typeI = function(t){
  k=seq(0,N)
  error=p[1:(N+1)]%*%(1-pnorm(t+k*mu*(shift)))
  return(error)
}
typeII = function(t){
  k=seq(0,N)
  error = p[1:(N+1)]%*%pnorm(t-k*mu*(1-shift))
  return(error)
}
Ivec2 = sapply(T,typeI)
IIvec2 = sapply(T,typeII)
mix_delta_new.df <- data.frame(Ivec1=Ivec2, IIvec1=IIvec2, line_type=line_type_val, alpha_val=rep(alpha_val, length(IIvec1)), data_pair=rep('5', length(IIvec1)))
mix_delta.df <- rbind(mix_delta.df, mix_delta_new.df)


# -1
shift = 1.0
typeI = function(t){
  k=seq(0,N)
  error=p[1:(N+1)]%*%(1-pnorm(t+k*mu*(shift)))
  return(error)
}
typeII = function(t){
  k=seq(0,N)
  error = p[1:(N+1)]%*%pnorm(t-k*mu*(1-shift))
  return(error)
}
Ivec2 = sapply(T,typeI)
IIvec2 = sapply(T,typeII)
mix_delta_new.df <- data.frame(Ivec1=Ivec2, IIvec1=IIvec2, line_type=line_type_val, alpha_val=rep(alpha_val, length(IIvec1)), data_pair=rep('6', length(IIvec1)))
mix_delta.df <- rbind(mix_delta.df, mix_delta_new.df)


# -1
shift = 1.0
typeI = function(t){
  k=seq(0,N)
  error=(1-pnorm(t+1*mu*(shift)))
  return(error)
}
typeII = function(t){
  k=seq(0,N)
  error = pnorm(t-1*mu*(1-shift))
  return(error)
}
Ivec2 = sapply(T,typeI)
IIvec2 = sapply(T,typeII)
mix_delta_new.df <- data.frame(Ivec1=Ivec2, IIvec1=IIvec2, line_type=2, alpha_val=rep(1, length(IIvec1)), data_pair=rep('7', length(IIvec1)))
mix_delta_gdp1.df <- rbind(mix_delta.df, mix_delta_new.df)



total_len <- 10001
half_len <- as.integer(total_len/2)
x_seq2 <- seq(-50,50,length=total_len)
mu_i <- 1
typeI_lowerbound <- 0 * x_seq2
typeII_lowerbound <- 0 * x_seq2
for (mu_i in seq(100)){
  typeI_lowerbound = typeI_lowerbound + (1-pnorm(x_seq2/2/mu_i+mu_i/2, mean=0)) * dbinom(mu_i,100,1/100)
  typeII_lowerbound = typeII_lowerbound + (pnorm(x_seq2/2/mu_i-mu_i/2, mean=0)) * dbinom(mu_i,100,1/100)
}
dbinom0 <- dbinom(0,100,1/100)
typeI_lowerbound <- typeI_lowerbound / (1-dbinom0)
typeII_lowerbound <- typeII_lowerbound / (1-dbinom0)
typeI_lowerbound[1:(half_len)] <- (1-dbinom0)*typeI_lowerbound[1:(half_len)] + 
  dbinom0 * (1-typeII_lowerbound[1:(half_len)])
typeII_lowerbound[(half_len+2):total_len] <- (1-dbinom0)*typeII_lowerbound[(half_len+2):total_len] + dbinom0 * (1-typeI_lowerbound[(half_len+2):total_len])
typeI_lowerbound <- c(0,typeI_lowerbound[1:(half_len)], typeI_lowerbound[(half_len+2):total_len])
typeII_lowerbound <- c(1,typeII_lowerbound[1:(half_len)], typeII_lowerbound[(half_len+2):total_len])

mix_delta_new.df <- data.frame(Ivec1=typeI_lowerbound, IIvec1=typeII_lowerbound, 
                               line_type=1,
                               alpha_val=rep(1, (total_len)), 
                               data_pair=rep('8', (total_len)))
mix_delta_lowerbound.df <- rbind(mix_delta_gdp1.df, mix_delta_new.df)

mix_delta_new2.df <- data.frame(Ivec1=c(1,1), IIvec1=c(1,1),
                                line_type=1,
                                alpha_val=rep(0, (2)),
                                data_pair=rep('9', (2)))
mix_delta_lowerbound.df <- rbind(mix_delta_gdp1.df, mix_delta_new.df, mix_delta_new2.df)


# The palette with black:
cbPalette <- c("#0000FF", "#E69F00", "#56B4E9", "#009E73", "#00FF00", "#0072B2", "#D55E00", "#000000", "#FF0000", "#CC79A7")
# cbPalette <- c("#AC8254", "#5AA800", "#0098E9", "#F29318", "#00FF00", "#0000FF", "#FF0000", "#000000")

alpha_val2 <- 0.36
override.alpha <- c(alpha_val2, alpha_val2, alpha_val2, alpha_val2, alpha_val2, alpha_val2, alpha_val2, 1, 1, 1)
override.linetype <- c(line_type_val,line_type_val,line_type_val,4,line_type_val,line_type_val,line_type_val,2,1,3)
p1 <- ggplot() + xlim(0,1) + xlab("Type I error")+ ylab("Type II error") +
  geom_line(data=mix_delta_lowerbound.df, aes(x=Ivec1, y=IIvec1, linetype=(line_type), alpha=factor(alpha_val), color=data_pair)) +
  scale_linetype_identity() +
  theme_classic() + scale_colour_manual(values=cbPalette, labels = c("0" = "a=0.0",
                                                                     "1" = "a=0.2",
                                                                     "2" = "a=0.4",
                                                                     "3" = "a=0.5",
                                                                     "4" = "a=0.6",
                                                                     "5" = "a=0.8",
                                                                     "6" = "a=1.0",
                                                                     "7" = "1-GDP",
                                                                     "8" = "Our lower bound",
                                                                     "9" = "1.125-GDP")) +
  scale_alpha_manual(values = c(0, 0.3, 1)) +
  theme(
    legend.position="right",
    legend.background = element_rect(fill='transparent'),
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
  ) + 
  guides(color=guide_legend(title="Tradeoff functions", ncol=2, bycol=TRUE,
                            override.aes = list(alpha = override.alpha, linetype = override.linetype)), 
         linetype="none", alpha="none") 








total_len <- 1000


x_seq <- seq(-3,3,length=total_len)
typeI_seq <- pnorm(-x_seq, mean=0.5)
typeII_seq <- pnorm(x_seq, mean=0.5)
typeI_seq2 <- pnorm(-x_seq, mean=0.5625)
typeII_seq2 <- pnorm(x_seq, mean=0.5625)
tradeoff.df <- data.frame(typeI_seq=c(typeI_seq, typeI_seq2), 
                          typeII_seq=c(typeII_seq, typeII_seq2),
                          privacy=rep(c("1-GDP", "1.125-GDP"), each=length(typeI_seq)),
                          lt=rep(c(2,3), each=length(typeI_seq)))
override.alpha <- c(1,1)
override.linetype <- c(2,3)
p2 <- ggplot(data=tradeoff.df, aes(x=typeI_seq, y=typeII_seq, group=privacy, color=privacy)) + 
  geom_line(aes(linetype=lt)) +
  scale_linetype_identity() +
  scale_color_manual(values=c("#000000", "#CC79A7")) +
  theme_classic() + theme(legend.position="right", panel.background = element_rect(fill='transparent'),
                          legend.background = element_rect(fill='transparent'),
                          plot.background = element_rect(fill='transparent', color=NA)) + 
  xlab('Type I error') + ylab('Type II error') + 
  guides(color=guide_legend(title="Tradeoff\nfunctions", ncol=1, 
                            override.aes = list(alpha = override.alpha, linetype = override.linetype)), 
         linetype="none", alpha="none") 





library(ggpubr)
ggarrange(p1, p2, widths = c(1,1), ncol=2, nrow=1, common.legend = TRUE, legend="right", labels="auto")
figure_width <- 7.5
figure_height <- 2.3
ggsave(paste("Figure3.pdf", sep=""), width=figure_width, height=figure_height)














