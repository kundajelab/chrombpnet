rm(list=ls())
library(ggplot2)
library(ggridges)
library(coin)

tf='spi1'
data=read.table(paste("5e-5.gm12878.",tf,".tsv",sep=''),header=TRUE,sep=',')
data$sig=data$pvalue<5e-5

data$Abs_Delta_Counts=abs(data$Alt_Minus_Ref)
data$Label=paste("Log Counts Delta")
# ggplot(data=data,
#        aes(x=data$Abs_Delta_Counts,
#            y=data$Label,
#            color=data$sig,
#            point_color=data$sig,
#            fill=data$sig))+
#   geom_density_ridges(
#     jittered_points=TRUE,
#     alpha=0.3,
#     scale=0.95,
#     rel_min_height=0.01,
#     point_shape = "|", point_size = 3, size = 0.25,
#     position = position_points_jitter(height = 0)) +
#   scale_y_discrete(expand = c(0, 0), name="") +
#   scale_x_continuous(expand = c(0, 0), name = "Abs(Log(Count(Alt))-Log(Count(Ref)))") +
#   scale_fill_manual(name="Sig bQTL",values = c('#e41a1ca0','#377eb8a0'), labels = c("FALSE", "TRUE")) +
#   scale_color_manual(values = c('#e41a1ca0','#377eb8a0'), guide = "none") +
#   scale_discrete_manual("point_color", values = c('#e41a1ca0','#377eb8a0'), guide = "none") +
#   coord_cartesian(clip = "off") +
#   guides(fill = guide_legend(
#     override.aes = list(
#       fill = c('#e41a1ca0','#377eb8a0'),
#       color = NA, point_color = NA)
#   )
#   ) +
#   ggtitle("Chrom-BPNET Count Preds for SPI1 bQTLs in GM12878") +
#   theme_ridges(center = TRUE)+
#   theme(legend.position = "top") 
data$sig=factor(data$sig)
sig_delta=data$Abs_Delta_Counts[data$sig==TRUE]
num_sig=length(sig_delta)
non_sig_delta=data$Abs_Delta_Counts[data$sig==FALSE]
num_non_sig=length(non_sig_delta)
#ks.test(sig_delta,non_sig_delta)
wilcox_test(Abs_Delta_Counts~sig,data=data,conf.level=0.99,conf.int=TRUE)

ggplot(data, aes(Abs_Delta_Counts, colour =sig)) +
  stat_ecdf(size=1)+
  scale_color_manual(values=c("#e41a1c","#377eb8a0"))+
  xlab("Abs(Log(Count(Alt))-Log(Count(Ref)))")+
  ylab("CDF")+
  ggtitle(paste(tf, ' p-val<5e-5 (n=',num_sig,' \n vs p-val>5e-5 (n=',num_non_sig,sep=''))+
  theme_bw(15)+
  theme(legend.position="top")
