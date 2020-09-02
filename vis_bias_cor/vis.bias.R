rm(list=ls())
library(ggplot2)
library(ggridges)
source("~/helpers.R")
data=read.table("bias.metrics.txt",header=TRUE,sep='\t')
dnase=data[data$Assay=="DNASE",]
atac=data[data$Assay=="ATAC",]

make_bias_dnase_plot <-function(data,metric){
  p=ggplot(data,
           aes(y=data$Model,
               x=data[[metric]],
               color=data$TestSet,
               point_color=data$TestSet,
               fill=data$TestSet))+
    geom_density_ridges(
      jittered_points=TRUE,
      scale=0.95,
      rel_min_height=0.01,
      point_shape = "|", point_size = 3, size = 0.25,
      position = position_points_jitter(height = 0)
    ) +
    scale_y_discrete(expand = c(0, 0), name="Model") +
    scale_x_continuous(expand = c(0, 0), name = metric) +
    scale_fill_manual(name="Smoothing",values = c("#D55E0050", "#FF000050","#0072B250"), labels = c("GM12878 Peaks", "K562 Peaks","Strided Genome")) +
    scale_color_manual(values = c("#D55E00", "#FF0000","#0072B2"), guide = "none") +
    scale_discrete_manual("point_color", values = c("#D55E00", "#FF0000", "#0072B2"), guide = "none") +
    coord_cartesian(clip = "off") +
    guides(fill = guide_legend(
      override.aes = list(
        fill = c("#D55E00A0","#FF0000A0", "#0072B2A0"),
        color = NA, point_color = NA)
    )
    ) +
    ggtitle("Bias, DNASE") +
    theme_ridges(center = TRUE)+
    theme(legend.position = "top") 
  return(p)
}

make_bias_atac_plot <-function(data,metric){
  p=ggplot(data,
           aes(y=data$Model,
               x=data[[metric]],
               color=data$TestSet,
               point_color=data$TestSet,
               fill=data$TestSet))+
    geom_density_ridges(
      jittered_points=TRUE,
      scale=0.95,
      rel_min_height=0.01,
      point_shape = "|", point_size = 3, size = 0.25,
      position = position_points_jitter(height = 0)
    ) +
    scale_y_discrete(expand = c(0, 0), name="Model") +
    scale_x_continuous(expand = c(0, 0), name = metric) +
    scale_fill_manual(name="Smoothing",values = c("#D55E0050", "#0072B250"), labels = c("GM12878 Peaks","Strided Genome")) +
    scale_color_manual(values = c("#D55E00","#0072B2"), guide = "none") +
    scale_discrete_manual("point_color", values = c("#D55E00", "#0072B2"), guide = "none") +
    coord_cartesian(clip = "off") +
    guides(fill = guide_legend(
      override.aes = list(
        fill = c("#D55E00A0", "#0072B2A0"),
        color = NA, point_color = NA)
    )
    ) +
    ggtitle("Bias ATAC") +
    theme_ridges(center = TRUE)+
    theme(legend.position = "top")
  return(p)
}


p1=make_bias_dnase_plot(dnase,"Pearson")
p2=make_bias_atac_plot(atac,"Pearson")

p3=make_bias_dnase_plot(dnase,"Spearman")
p4=make_bias_atac_plot(atac,"Spearman")

p5=make_bias_dnase_plot(dnase,"MSE")
p6=make_bias_atac_plot(atac,"MSE")

p7=make_bias_dnase_plot(dnase,"MeanJSD")
p8=make_bias_atac_plot(atac,"MeanJSD")

p9=make_bias_dnase_plot(dnase,"MeanMNNLL")
p10=make_bias_atac_plot(atac,"MeanMNNLL")

#remove mse outliers 
p11=make_bias_dnase_plot(dnase,"MSE")+xlim(c(-1,2))
p12=make_bias_atac_plot(atac,"MSE")+xlim(c(-1,1))