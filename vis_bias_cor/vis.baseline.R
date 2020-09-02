rm(list=ls())
library(ggplot2)
library(ggridges)
source("~/helpers.R")
data=read.table("baseline.metrics.txt",header=TRUE,sep='\t')
atac=data[data$Assay=="ATAC",]
dnase=data[data$Assay=="DNASE",]

make_baseline_plot <-function(data,metric,cur_title){
p=ggplot(data,
       aes(y=data$CellLine,
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
  scale_y_discrete(expand = c(0, 0), name="Cell Line") +
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
    
p1=make_baseline_plot(dnase,'Pearson','DNASE, Baseline')
p2=make_baseline_plot(dnase,'Spearman','DNASE, Baseline')
p3=make_baseline_plot(dnase,'MSE','DNASE, Baseline')
p4=make_baseline_plot(dnase,'MeanJSD','DNASE, Baseline')
p5=make_baseline_plot(dnase,'MeanMNNLL','DNASE, Baseline')

p6=make_baseline_plot(atac,'Pearson','ATAC, Baseline')
p7=make_baseline_plot(atac,'Spearman','ATAC, Baseline')
p8=make_baseline_plot(atac,'MSE','ATAC, Baseline')
p9=make_baseline_plot(atac,'MeanJSD','ATAC, Baseline')
p10=make_baseline_plot(atac,'MeanMNNLL','ATAC, Baseline')

multiplot(p1,p6,p2,p7,p3,p8,p4,p9,p5,p10,cols=5)
