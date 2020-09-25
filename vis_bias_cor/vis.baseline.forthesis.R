rm(list=ls())
library(reshape2)
library(ggplot2)
library(ggridges)
source("~/helpers.R")
data=read.table("baseline.metrics.txt",header=TRUE,sep='\t')
data$CellLine=factor(data$CellLine,levels=c("K562","GM12878","IMR90","HEPG2","H1ESC"))
data=data[data$Smoothing=="Both",]
atac=data[data$Assay=="ATAC",]
dnase=data[data$Assay=="DNASE",]
p1=ggplot(data=dnase,
         aes(x=CellLine,
             y=Spearman,
             group=CellLine,
             fill=CellLine))+
  geom_boxplot()+
  geom_point(position=position_jitterdodge())+
  scale_fill_manual(values=c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00'))+
  ylim(0.4,0.9)+
  xlab("Cell Line")+
  ylab("Spearman R, Test Set")+
  theme_bw(20)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


p2=ggplot(data=atac,
          aes(x=CellLine,
              y=Spearman,
              group=CellLine,
              fill=CellLine))+
  geom_boxplot()+
  geom_point(position=position_jitterdodge())+
  scale_fill_manual(values=c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00'))+
  ylim(0.4,0.9)+
  xlab("Cell Line")+
  ylab("Spearman R, Test Set")+
  theme_bw(20)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


jsd_dnase=melt(read.table("dnase.smoothed.both.jsd.all.txt",header=TRUE,sep='\t'))
jsd_dnase$variable=factor(jsd_dnase$variable,levels=c("K562","GM12878","IMR90","HEPG2","H1ESC"))

jsd_atac=melt(read.table("atac.smoothed.both.jsd.all.txt",header=TRUE,sep='\t'))
jsd_atac$variable=factor(jsd_atac$variable,levels=c("GM12878","IMR90","HEPG2","H1ESC"))

#mnll_dnase=melt(read.table("dnase.smoothed.both.jsd.all.txt",header=TRUE,sep='\t'))
#mnll_atac=melt(read.table("atac.smoothed.both.mnll.all.txt",header=TRUE,sep='\t'))

make_baseline_plot <-function(data,cur_title){
p=ggplot(data,
       aes(y=data$variable,
           x=data$value,
           fill=data$variable))+
  geom_density_ridges(
    jittered_points=TRUE,
    scale=0.95,
    rel_min_height=0.01,
    point_shape = "|", point_size = 1, size = 0.1,
    position = position_points_jitter(height = 0)
  ) +
  scale_y_discrete(expand = c(0, 0), name="Cell Line") +
  scale_x_continuous(expand = c(0, 0), name = "JSD") +
  scale_fill_manual(name="Cell Line",values = c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00')) +
  coord_cartesian(clip = "off") +
  ggtitle(cur_title) +
  theme_ridges(center = TRUE) +
  theme(legend.position = "top")
return(p)
}
    
p3=make_baseline_plot(jsd_dnase,'JSD DNASE')
p4=make_baseline_plot(jsd_atac,'JSD ATAC')
