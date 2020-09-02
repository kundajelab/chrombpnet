rm(list=ls())
library(ggplot2)
library(ggridges)
source("~/helpers.R")
data=read.table("bias.corrected.metrics.txt",header=TRUE,sep='\t')
make_bias_corrected_plot <-function(data,metric,cur_title){
  p=ggplot(data,
           aes(y=data$Model,
               x=data[[metric]],
               color=data$Smoothing,
               point_color=data$Smoothing,
               fill=data$Smoothing))+
    geom_density_ridges(
      jittered_points=TRUE,
      scale=0.95,
      rel_min_height=0.01,
      point_shape = "|", point_size = 3, size = 0.25,
      position = position_points_jitter(height = 0)
    ) +
    scale_y_discrete(expand = c(0, 0), name="Model") +
    scale_x_continuous(expand = c(0, 0), name = metric) +
    scale_fill_manual(name="Smoothing",values = c("#D55E0050", "#FF000050","#0072B250"), labels = c("Both", "Labels","None")) +
    scale_color_manual(values = c("#D55E00", "#FF0000","#0072B2"), guide = "none") +
    scale_discrete_manual("point_color", values = c("#D55E00", "#FF0000", "#0072B2"), guide = "none") +
    coord_cartesian(clip = "off") +
    guides(fill = guide_legend(
      override.aes = list(
        fill = c("#D55E00A0","#FF0000A0", "#0072B2A0"),
        color = NA, point_color = NA)
    )
    ) +
    ggtitle(cur_title) +
    theme_ridges(center = TRUE)+
    theme(legend.position = "top") 
  return(p)
}

#cell="GM12878"
cell="K562"
assay="DNASE"
p1=make_bias_corrected_plot(data[data$CellLine==cell & data$Assay==assay,], "Pearson",paste('Pearson',cell,assay,"Corrected"))
p2=make_bias_corrected_plot(data[data$CellLine==cell & data$Assay==assay,],"Spearman",paste('Spearman',cell,assay,"Corrected"))
p3=make_bias_corrected_plot(data[data$CellLine==cell & data$Assay==assay,],"MSE",paste('MSE',cell,assay,"Corrected"))
p4=make_bias_corrected_plot(data[data$CellLine==cell & data$Assay==assay,],"MeanJSD",paste('MeanJSD',cell,assay,"Corrected"))
p5=make_bias_corrected_plot(data[data$CellLine==cell & data$Assay==assay,],"MeanMNNLL",paste('MeanMNNLL',cell,assay,"Corrected"))
multiplot(p1,p4,p2,p5,p3,cols=3)

