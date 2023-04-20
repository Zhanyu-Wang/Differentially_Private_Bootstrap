library(ggplot2)
library(ggpubr)
cbPalette <- c("#999999","#E69F00","#56B4E9","#009E73",
               "#F0E442","#0072B2","#D55E00","#CC79A7")
n <- 100
clean_x <- seq(0,4,1)
noisy_x <- seq(-3,4.5,length=1000)
clean_hist <- dbinom(clean_x, n, 1/n)
component_distI <- dnorm(noisy_x, 0, 1)
component_distII_0 <- dnorm(noisy_x, 0, 1) * dbinom(0, n, 1/n)
component_distII_1 <- dnorm(noisy_x, 1, 1) * dbinom(1, n, 1/n) 
component_distII_2 <- dnorm(noisy_x, 2, 1) * dbinom(2, n, 1/n) 
component_distII_3 <- dnorm(noisy_x, 3, 1) * dbinom(3, n, 1/n) 
component_distII_4 <- dnorm(noisy_x, 4, 1) * dbinom(4, n, 1/n) 
component_distII <- component_distII_0 + component_distII_1 + component_distII_2 + component_distII_3 + component_distII_4
hist.df <- data.frame(clean_x, clean_hist)
density0.df <- data.frame(noisy_x=rep(noisy_x, 2), 
                         density=c(component_distI, component_distII),
                         comp_idx=rep(seq(1,2),each=1000))
density.df <- data.frame(noisy_x=rep(noisy_x, 5), 
                         density=c(component_distII_0, 
                         component_distII_1, component_distII_2, 
                         component_distII_3, component_distII_4) / dnorm(0),
                         comp_idx=rep(seq(1,5),each=1000))


total_len <- 10001

x_seq <- seq(-3,3,length=total_len)
typeI_seq0 <- 1-pnorm(x_seq, mean=0)
typeII_seq0 <- pnorm(x_seq, mean=0)
typeI_seq1 <- 1-pnorm(x_seq, mean=0)
typeII_seq1 <- pnorm(x_seq, mean=1)
typeI_seq2 <- 1-pnorm(x_seq, mean=0)
typeII_seq2 <- pnorm(x_seq, mean=2)
typeI_seq3 <- 1-pnorm(x_seq, mean=0)
typeII_seq3 <- pnorm(x_seq, mean=3)
tradeoff.df <- data.frame(typeI_seq=c(0,typeI_seq1,0,typeI_seq2,0,typeI_seq3), 
                          typeII_seq=c(1,typeII_seq1,1,typeII_seq2,1,typeII_seq3),
                          privacy=c(rep(c("1-GDP","2-GDP","3-GDP"), each=1+length(typeI_seq0))))

rej_x <- c(-1)
rej_typeI <- pnorm(-rej_x, mean=0.5)
rej_typeII <- pnorm(rej_x, mean=0.5)
Reject_H0 <- c('2')
rej_rules.df <- data.frame(rej_x=rej_x, rej_typeI=rej_typeI, rej_typeII=rej_typeII, Reject_H0=Reject_H0)

# only for lower bound
half_len <- as.integer(total_len/2)
x_seq2 <- seq(-50,50,length=total_len)
mu_i <- 1
typeI_lowerbound <- 0 * x_seq2
typeII_lowerbound <- 0 * x_seq2
for (mu_i in c(1,2,3)){
typeI_lowerbound = typeI_lowerbound + (1-pnorm(x_seq2/2/mu_i+mu_i/2, mean=0)) * 1/3# (dbinom(mu_i,100,1/100) / (dbinom(1,100,1/100)+dbinom(2,100,1/100)+dbinom(3,100,1/100)))
typeII_lowerbound = typeII_lowerbound + (pnorm(x_seq2/2/mu_i-mu_i/2, mean=0)) * 1/3# (dbinom(mu_i,100,1/100) / (dbinom(1,100,1/100)+dbinom(2,100,1/100)+dbinom(3,100,1/100)))
}
typeI_lowerbound <- c(0,typeI_lowerbound[1:(half_len)], typeI_lowerbound[(half_len+2):total_len])
typeII_lowerbound <- c(1,typeII_lowerbound[1:(half_len)], typeII_lowerbound[(half_len+2):total_len])


mix_delta_lowerbound.df <- data.frame(Ivec1=typeI_lowerbound, IIvec1=typeII_lowerbound, 
                               alpha_val=rep(1, (total_len)), 
                               data_pair=rep('7', (total_len)))

rej_x = 3
typeI_lowerbound <- 0 
typeII_lowerbound <- 0 
for (mu_i in c(1,2,3)){
  typeI_lowerbound = typeI_lowerbound + (1-pnorm(rej_x/2/mu_i+mu_i/2, mean=0)) * 1/3# (dbinom(mu_i,100,1/100) / (dbinom(1,100,1/100)+dbinom(2,100,1/100)+dbinom(3,100,1/100)))
  typeII_lowerbound = typeII_lowerbound + (pnorm(rej_x/2/mu_i-mu_i/2, mean=0)) * 1/3# (dbinom(mu_i,100,1/100) / (dbinom(1,100,1/100)+dbinom(2,100,1/100)+dbinom(3,100,1/100)))
}

mix_delta_lowerbound_point.df <- data.frame(Ivec1=typeI_lowerbound, IIvec1=typeII_lowerbound, 
                                            alpha_val=1, 
                                            data_pair="7")

lb3 <- ggplot() + xlim(0,1) + ylim(0,1) + xlab("Type I error")+ ylab("Type II error") +
  geom_line(data=mix_delta_lowerbound.df, aes(x=Ivec1, y=IIvec1, color=data_pair)) +
  geom_point(data=mix_delta_lowerbound_point.df, aes(x=Ivec1, y=IIvec1, color=data_pair), size=2, pch=21, fill='white') +
  theme_classic() + 
  theme(
    legend.position="none",
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
  )

rej_x = 2
typeI_lowerbound <- 0 
typeII_lowerbound <- 0 
for (mu_i in c(1,2,3)){
  typeI_lowerbound = typeI_lowerbound + (1-pnorm(rej_x/2/mu_i+mu_i/2, mean=0)) * 1/3# (dbinom(mu_i,100,1/100) / (dbinom(1,100,1/100)+dbinom(2,100,1/100)+dbinom(3,100,1/100)))
  typeII_lowerbound = typeII_lowerbound + (pnorm(rej_x/2/mu_i-mu_i/2, mean=0)) * 1/3# (dbinom(mu_i,100,1/100) / (dbinom(1,100,1/100)+dbinom(2,100,1/100)+dbinom(3,100,1/100)))
}

mix_delta_lowerbound_point.df <- data.frame(Ivec1=typeI_lowerbound, IIvec1=typeII_lowerbound, 
                                            alpha_val=1, 
                                            data_pair="7")

lb2 <- ggplot() + xlim(0,1) + ylim(0,1) + xlab("Type I error")+ ylab("Type II error") +
  geom_line(data=mix_delta_lowerbound.df, aes(x=Ivec1, y=IIvec1, color=data_pair)) +
  geom_point(data=mix_delta_lowerbound_point.df, aes(x=Ivec1, y=IIvec1, color=data_pair), size=2, pch=21, fill='white') +
  theme_classic() + 
  theme(
    legend.position="none",
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
  )

rej_x = 1
typeI_lowerbound <- 0 
typeII_lowerbound <- 0 
for (mu_i in c(1,2,3)){
  typeI_lowerbound = typeI_lowerbound + (1-pnorm(rej_x/2/mu_i+mu_i/2, mean=0)) * 1/3# (dbinom(mu_i,100,1/100) / (dbinom(1,100,1/100)+dbinom(2,100,1/100)+dbinom(3,100,1/100)))
  typeII_lowerbound = typeII_lowerbound + (pnorm(rej_x/2/mu_i-mu_i/2, mean=0)) * 1/3# (dbinom(mu_i,100,1/100) / (dbinom(1,100,1/100)+dbinom(2,100,1/100)+dbinom(3,100,1/100)))
}

mix_delta_lowerbound_point.df <- data.frame(Ivec1=typeI_lowerbound, IIvec1=typeII_lowerbound, 
                                            alpha_val=1, 
                                            data_pair="7")

lb1 <- ggplot() + xlim(0,1) + ylim(0,1) + xlab("Type I error")+ ylab("Type II error") +
  geom_line(data=mix_delta_lowerbound.df, aes(x=Ivec1, y=IIvec1, color=data_pair)) +
  geom_point(data=mix_delta_lowerbound_point.df, aes(x=Ivec1, y=IIvec1, color=data_pair), size=2, pch=21, fill='white') +
  theme_classic() + 
  theme(
    legend.position="none",
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
  )

rej_x = 0
typeI_lowerbound <- 0 
typeII_lowerbound <- 0 
for (mu_i in c(1,2,3)){
  typeI_lowerbound = typeI_lowerbound + (1-pnorm(rej_x/2/mu_i+mu_i/2, mean=0)) * 1/3# (dbinom(mu_i,100,1/100) / (dbinom(1,100,1/100)+dbinom(2,100,1/100)+dbinom(3,100,1/100)))
  typeII_lowerbound = typeII_lowerbound + (pnorm(rej_x/2/mu_i-mu_i/2, mean=0)) * 1/3# (dbinom(mu_i,100,1/100) / (dbinom(1,100,1/100)+dbinom(2,100,1/100)+dbinom(3,100,1/100)))
}

mix_delta_lowerbound_point.df <- data.frame(Ivec1=typeI_lowerbound, IIvec1=typeII_lowerbound, 
                                            alpha_val=1, 
                                            data_pair="7")

lb0 <- ggplot() + xlim(0,1) + ylim(0,1) + xlab("Type I error")+ ylab("Type II error") +
  geom_line(data=mix_delta_lowerbound.df, aes(x=Ivec1, y=IIvec1, color=data_pair)) +
  geom_point(data=mix_delta_lowerbound_point.df, aes(x=Ivec1, y=IIvec1, color=data_pair), size=2, pch=21, fill='white') +
  theme_classic() + 
  theme(
    legend.position="none",
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
  )

x_seq <- seq(-3,3,length=total_len)
typeI_seq0 <- 1-pnorm(x_seq, mean=0)
typeII_seq0 <- pnorm(x_seq, mean=0)
typeI_seq1 <- 1-pnorm(x_seq, mean=0)
typeII_seq1 <- pnorm(x_seq, mean=1)
typeI_seq2 <- 1-pnorm(x_seq, mean=0)
typeII_seq2 <- pnorm(x_seq, mean=2)
typeI_seq3 <- 1-pnorm(x_seq, mean=0)
typeII_seq3 <- pnorm(x_seq, mean=3)
tradeoff.df <- data.frame(typeI_seq=c(0,typeI_seq1,0,typeI_seq2,0,typeI_seq3), 
                          typeII_seq=c(1,typeII_seq1,1,typeII_seq2,1,typeII_seq3),
                          privacy=c(rep(c("1-GDP","2-GDP","3-GDP"), each=1+length(typeI_seq0))))


rej_x = 3
rej_typeI_seq0 <- 1-pnorm(rej_x/2/mu_i+mu_i/2, mean=0)
rej_typeII_seq0 <- pnorm(rej_x/2/mu_i-mu_i/2, mean=0)

tangent_param <- seq(-1, 1, length=100)
mu_i <- 1
rej_typeI_seq1 <- 1-pnorm(rej_x/2/mu_i+mu_i/2, mean=0)
rej_typeII_seq1 <- pnorm(rej_x/2/mu_i-mu_i/2, mean=0)
rej_typeI_seq1_x <- rej_typeI_seq1 + tangent_param * dnorm(rej_x/2/mu_i+mu_i/2, mean=0)
rej_typeII_seq1_x <- rej_typeII_seq1 - tangent_param * dnorm(rej_x/2/mu_i-mu_i/2, mean=0)
mu_i <- 2
rej_typeI_seq2 <- 1-pnorm(rej_x/2/mu_i+mu_i/2, mean=0)
rej_typeII_seq2 <- pnorm(rej_x/2/mu_i-mu_i/2, mean=0)
rej_typeI_seq2_x <- rej_typeI_seq2 + tangent_param * dnorm(rej_x/2/mu_i+mu_i/2, mean=0)
rej_typeII_seq2_x <- rej_typeII_seq2 - tangent_param * dnorm(rej_x/2/mu_i-mu_i/2, mean=0)
mu_i <- 3
rej_typeI_seq3 <- 1-pnorm(rej_x/2/mu_i+mu_i/2, mean=0)
rej_typeII_seq3 <- pnorm(rej_x/2/mu_i-mu_i/2, mean=0)
rej_typeI_seq3_x <- rej_typeI_seq3 + tangent_param * dnorm(rej_x/2/mu_i+mu_i/2, mean=0)
rej_typeII_seq3_x <- rej_typeII_seq3 - tangent_param * dnorm(rej_x/2/mu_i-mu_i/2, mean=0)
rej_typeI <- c(rej_typeI_seq1, rej_typeI_seq2, rej_typeI_seq3)
rej_typeII <- c(rej_typeII_seq1, rej_typeII_seq2, rej_typeII_seq3)
rej_typeI_x <- c(rej_typeI_seq1_x, rej_typeI_seq2_x, rej_typeI_seq3_x)
rej_typeII_x <- c(rej_typeII_seq1_x, rej_typeII_seq2_x, rej_typeII_seq3_x)

Reject_H0 <- c("1-GDP","2-GDP","3-GDP")
Reject_H0_x <- rep(c("1-GDP","2-GDP","3-GDP"), each=100)
rej_rules.df3 <- data.frame(rej_x=rej_x, rej_typeI=rej_typeI, rej_typeII=rej_typeII, privacy=Reject_H0)
tangent_line.df3 <- data.frame(rej_typeI=rej_typeI_x, rej_typeII=rej_typeII_x, privacy=Reject_H0_x)


mix3 <- ggplot(data=tradeoff.df, aes(x=typeI_seq, y=typeII_seq, group=privacy, color=privacy)) + 
  geom_line()  + xlim(0,1) + ylim(0,1) +
  geom_line(data=tangent_line.df3, aes(x=rej_typeI, y=rej_typeII, group=privacy, color=privacy), linetype = "dashed") +
  geom_point(data=rej_rules.df3, aes(x=rej_typeI, y=rej_typeII, color=privacy), size=2, pch=21, fill='white') +
  theme_classic() + scale_colour_manual(values=cbPalette) +
  theme(legend.position="top", panel.background = element_rect(fill='transparent'),
        legend.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA)) + 
  xlab('Type I Error') + ylab('Type II Error') +
  guides(color=guide_legend(nrow=1,byrow=TRUE))


rej_x = 2
rej_typeI_seq0 <- 1-pnorm(rej_x/2/mu_i+mu_i/2, mean=0)
rej_typeII_seq0 <- pnorm(rej_x/2/mu_i-mu_i/2, mean=0)

tangent_param <- seq(-0.7, 0.7, length=100)
mu_i <- 1
rej_typeI_seq1 <- 1-pnorm(rej_x/2/mu_i+mu_i/2, mean=0)
rej_typeII_seq1 <- pnorm(rej_x/2/mu_i-mu_i/2, mean=0)
rej_typeI_seq1_x <- rej_typeI_seq1 + tangent_param * dnorm(rej_x/2/mu_i+mu_i/2, mean=0)
rej_typeII_seq1_x <- rej_typeII_seq1 - tangent_param * dnorm(rej_x/2/mu_i-mu_i/2, mean=0)
mu_i <- 2
rej_typeI_seq2 <- 1-pnorm(rej_x/2/mu_i+mu_i/2, mean=0)
rej_typeII_seq2 <- pnorm(rej_x/2/mu_i-mu_i/2, mean=0)
rej_typeI_seq2_x <- rej_typeI_seq2 + tangent_param * dnorm(rej_x/2/mu_i+mu_i/2, mean=0)
rej_typeII_seq2_x <- rej_typeII_seq2 - tangent_param * dnorm(rej_x/2/mu_i-mu_i/2, mean=0)
mu_i <- 3
rej_typeI_seq3 <- 1-pnorm(rej_x/2/mu_i+mu_i/2, mean=0)
rej_typeII_seq3 <- pnorm(rej_x/2/mu_i-mu_i/2, mean=0)
rej_typeI_seq3_x <- rej_typeI_seq3 + tangent_param * dnorm(rej_x/2/mu_i+mu_i/2, mean=0)
rej_typeII_seq3_x <- rej_typeII_seq3 - tangent_param * dnorm(rej_x/2/mu_i-mu_i/2, mean=0)
rej_typeI <- c(rej_typeI_seq1, rej_typeI_seq2, rej_typeI_seq3)
rej_typeII <- c(rej_typeII_seq1, rej_typeII_seq2, rej_typeII_seq3)
rej_typeI_x <- c(rej_typeI_seq1_x, rej_typeI_seq2_x, rej_typeI_seq3_x)
rej_typeII_x <- c(rej_typeII_seq1_x, rej_typeII_seq2_x, rej_typeII_seq3_x)

Reject_H0 <- c("1-GDP","2-GDP","3-GDP")
Reject_H0_x <- rep(c("1-GDP","2-GDP","3-GDP"), each=100)
rej_rules.df2 <- data.frame(rej_x=rej_x, rej_typeI=rej_typeI, rej_typeII=rej_typeII, privacy=Reject_H0)
tangent_line.df2 <- data.frame(rej_typeI=rej_typeI_x, rej_typeII=rej_typeII_x, privacy=Reject_H0_x)


mix2 <- ggplot(data=tradeoff.df, aes(x=typeI_seq, y=typeII_seq, group=privacy, color=privacy)) + 
  geom_line()  + xlim(0,1) + ylim(0,1) +
  geom_line(data=tangent_line.df2, aes(x=rej_typeI, y=rej_typeII, group=privacy, color=privacy), linetype = "dashed") +
  geom_point(data=rej_rules.df2, aes(x=rej_typeI, y=rej_typeII, color=privacy), size=2, pch=21, fill='white') +
  theme_classic() + scale_colour_manual(values=cbPalette) +
  theme(legend.position="top", panel.background = element_rect(fill='transparent'),
        legend.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA)) + 
  xlab('Type I Error') + ylab('Type II Error') +
  guides(color=guide_legend(nrow=1,byrow=TRUE))


rej_x = 1
rej_typeI_seq0 <- 1-pnorm(rej_x/2/mu_i+mu_i/2, mean=0)
rej_typeII_seq0 <- pnorm(rej_x/2/mu_i-mu_i/2, mean=0)

tangent_param <- seq(-0.7, 0.7, length=100)
mu_i <- 1
rej_typeI_seq1 <- 1-pnorm(rej_x/2/mu_i+mu_i/2, mean=0)
rej_typeII_seq1 <- pnorm(rej_x/2/mu_i-mu_i/2, mean=0)
rej_typeI_seq1_x <- rej_typeI_seq1 + tangent_param * dnorm(rej_x/2/mu_i+mu_i/2, mean=0)
rej_typeII_seq1_x <- rej_typeII_seq1 - tangent_param * dnorm(rej_x/2/mu_i-mu_i/2, mean=0)
mu_i <- 2
rej_typeI_seq2 <- 1-pnorm(rej_x/2/mu_i+mu_i/2, mean=0)
rej_typeII_seq2 <- pnorm(rej_x/2/mu_i-mu_i/2, mean=0)
rej_typeI_seq2_x <- rej_typeI_seq2 + tangent_param * dnorm(rej_x/2/mu_i+mu_i/2, mean=0)
rej_typeII_seq2_x <- rej_typeII_seq2 - tangent_param * dnorm(rej_x/2/mu_i-mu_i/2, mean=0)
mu_i <- 3
rej_typeI_seq3 <- 1-pnorm(rej_x/2/mu_i+mu_i/2, mean=0)
rej_typeII_seq3 <- pnorm(rej_x/2/mu_i-mu_i/2, mean=0)
rej_typeI_seq3_x <- rej_typeI_seq3 + tangent_param * dnorm(rej_x/2/mu_i+mu_i/2, mean=0)
rej_typeII_seq3_x <- rej_typeII_seq3 - tangent_param * dnorm(rej_x/2/mu_i-mu_i/2, mean=0)
rej_typeI <- c(rej_typeI_seq1, rej_typeI_seq2, rej_typeI_seq3)
rej_typeII <- c(rej_typeII_seq1, rej_typeII_seq2, rej_typeII_seq3)
rej_typeI_x <- c(rej_typeI_seq1_x, rej_typeI_seq2_x, rej_typeI_seq3_x)
rej_typeII_x <- c(rej_typeII_seq1_x, rej_typeII_seq2_x, rej_typeII_seq3_x)

Reject_H0 <- c("1-GDP","2-GDP","3-GDP")
Reject_H0_x <- rep(c("1-GDP","2-GDP","3-GDP"), each=100)
rej_rules.df1 <- data.frame(rej_x=rej_x, rej_typeI=rej_typeI, rej_typeII=rej_typeII, privacy=Reject_H0)
tangent_line.df1 <- data.frame(rej_typeI=rej_typeI_x, rej_typeII=rej_typeII_x, privacy=Reject_H0_x)


mix1 <- ggplot(data=tradeoff.df, aes(x=typeI_seq, y=typeII_seq, group=privacy, color=privacy)) + 
  geom_line()  + xlim(0,1) + ylim(0,1) +
  geom_line(data=tangent_line.df1, aes(x=rej_typeI, y=rej_typeII, group=privacy, color=privacy), linetype = "dashed") +
  geom_point(data=rej_rules.df1, aes(x=rej_typeI, y=rej_typeII, color=privacy), size=2, pch=21, fill='white') +
  theme_classic() + scale_colour_manual(values=cbPalette) +
  theme(legend.position="top", panel.background = element_rect(fill='transparent'),
        legend.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA)) + 
  xlab('Type I Error') + ylab('Type II Error') +
  guides(color=guide_legend(nrow=1,byrow=TRUE))


rej_x = 0
rej_typeI_seq0 <- 1-pnorm(rej_x/2/mu_i+mu_i/2, mean=0)
rej_typeII_seq0 <- pnorm(rej_x/2/mu_i-mu_i/2, mean=0)

tangent_param <- seq(-0.7, 0.7, length=100)
mu_i <- 1
rej_typeI_seq1 <- 1-pnorm(rej_x/2/mu_i+mu_i/2, mean=0)
rej_typeII_seq1 <- pnorm(rej_x/2/mu_i-mu_i/2, mean=0)
rej_typeI_seq1_x <- rej_typeI_seq1 + tangent_param * dnorm(rej_x/2/mu_i+mu_i/2, mean=0)
rej_typeII_seq1_x <- rej_typeII_seq1 - tangent_param * dnorm(rej_x/2/mu_i-mu_i/2, mean=0)
mu_i <- 2
rej_typeI_seq2 <- 1-pnorm(rej_x/2/mu_i+mu_i/2, mean=0)
rej_typeII_seq2 <- pnorm(rej_x/2/mu_i-mu_i/2, mean=0)
rej_typeI_seq2_x <- rej_typeI_seq2 + tangent_param * dnorm(rej_x/2/mu_i+mu_i/2, mean=0)
rej_typeII_seq2_x <- rej_typeII_seq2 - tangent_param * dnorm(rej_x/2/mu_i-mu_i/2, mean=0)
mu_i <- 3
rej_typeI_seq3 <- 1-pnorm(rej_x/2/mu_i+mu_i/2, mean=0)
rej_typeII_seq3 <- pnorm(rej_x/2/mu_i-mu_i/2, mean=0)
rej_typeI_seq3_x <- rej_typeI_seq3 + tangent_param * dnorm(rej_x/2/mu_i+mu_i/2, mean=0)
rej_typeII_seq3_x <- rej_typeII_seq3 - tangent_param * dnorm(rej_x/2/mu_i-mu_i/2, mean=0)
rej_typeI <- c(rej_typeI_seq1, rej_typeI_seq2, rej_typeI_seq3)
rej_typeII <- c(rej_typeII_seq1, rej_typeII_seq2, rej_typeII_seq3)
rej_typeI_x <- c(rej_typeI_seq1_x, rej_typeI_seq2_x, rej_typeI_seq3_x)
rej_typeII_x <- c(rej_typeII_seq1_x, rej_typeII_seq2_x, rej_typeII_seq3_x)

Reject_H0 <- c("1-GDP","2-GDP","3-GDP")
Reject_H0_x <- rep(c("1-GDP","2-GDP","3-GDP"), each=100)
rej_rules.df <- data.frame(rej_x=rej_x, rej_typeI=rej_typeI, rej_typeII=rej_typeII, privacy=Reject_H0)
tangent_line.df <- data.frame(rej_typeI=rej_typeI_x, rej_typeII=rej_typeII_x, privacy=Reject_H0_x)


mix0 <- ggplot(data=tradeoff.df, aes(x=typeI_seq, y=typeII_seq, group=privacy, color=privacy)) + 
  geom_line()  + xlim(0,1) + ylim(0,1) +
  geom_line(data=tangent_line.df, aes(x=rej_typeI, y=rej_typeII, group=privacy, color=privacy), linetype = "dashed") +
  geom_point(data=rej_rules.df, aes(x=rej_typeI, y=rej_typeII, color=privacy), size=2, pch=21, fill='white') +
  theme_classic() + scale_colour_manual(values=cbPalette) +
  theme(legend.position="top", panel.background = element_rect(fill='transparent'),
        legend.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA)) + 
  xlab('Type I Error') + ylab('Type II Error') +
  guides(color=guide_legend(nrow=1,byrow=TRUE))


ggarrange(mix3, mix2, mix1, mix0, lb3, lb2, lb1, lb0, 
          ncol = 4, nrow = 2,  align = "hv", 
          widths = c(1, 1, 1, 1), heights = c(2, 2),
          common.legend = TRUE) %>%
  ggexport(filename = "Figure2.pdf", width=8, height=4)

