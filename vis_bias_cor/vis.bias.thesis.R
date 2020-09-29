rm(list=ls())
library(ggplot2)
library(ggridges)
source("~/helpers.R")
data=read.table("bias.metrics.thesis.txt",header=TRUE,sep='\t')
data$Model=factor(data$Model)
data$Assay=factor(data$Assay)
p1=ggplot(data=data,
       aes(x=Model, 
           y=Spearman,
           fill=Assay))+
  geom_boxplot(position='dodge')+
  geom_jitter(shape=21,position=position_jitterdodge())+
  scale_fill_manual(values=c("#e41a1c","#377eb8"))+
  theme_bw(20)+
  coord_flip()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(0,1)+
  ylab("Spearman R\nfor Counts Pred")+
  theme(legend.position="bottom")

p2=ggplot(data=data,
          aes(x=Model, 
              y=MeanJSD,
              fill=Assay))+
  geom_boxplot()+
  geom_jitter(shape=21,position=position_jitterdodge())+
  scale_fill_manual(values=c("#e41a1c","#377eb8"))+
  theme_bw(20)+
  coord_flip()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(0,1)+
  ylab("Mean Jensen-Shannon Distance\nfor Profile Pred")+
  theme(legend.position="bottom")

multiplot(p1,p2,cols=2)