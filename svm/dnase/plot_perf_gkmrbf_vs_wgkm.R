rm(list=ls())
library(ggplot2)
data=read.table("aggregate_performance_wgkm_vs_gkmrbf.tsv",header=TRUE,sep='\t')
data=data[data$Metric %in% c("auprc", 
                             "auroc", 
                             "balanced_accuracy", 
                             "negative_accuracy",
                             "positive_accuracy",
                             "recall_at_fdr_50",
                             "recall_at_fdr_20"),]
attach(data)
ggplot(data, aes(x = Fold,y=Val,group=Model,fill=Model)) +
  geom_bar(stat='identity',position='dodge') +
  facet_grid(rows=vars(Metric),cols=vars(Task))+
  scale_fill_manual(values=c("#e41a1c",'#377eb8'))+
  scale_x_continuous(name="Fold",breaks=seq(0,9),labels=seq(0,9))

roc=data[data$Metric %in% c('auroc'),]
wgkm=roc[roc$Model=="wgkm",]
wgkm$WGKM=wgkm$Val
wgkm$Metric=NULL
wgkm$Model=NULL
wgkm$Val=NULL

gkmrbf=roc[roc$Model=="gkmrbf",]
gkmrbf$GKMRBF=gkmrbf$Val
gkmrbf$Metric=NULL
gkmrbf$Model=NULL
gkmrbf$Val=NULL

merged=merge(gkmrbf,wgkm)
ggplot(merged,aes(x=merged$GKMRBF,y=merged$WGKM,color=merged$Task,label=merged$Fold))+
  geom_point()+
  geom_text()+
  scale_color_manual(values=c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00'))+
  geom_abline()+
  xlab("GKMRBF")+
  ylab("WGKM")+
  ggtitle('AUROC')+
  xlim(0.7,0.9)+
  ylim(0.7,0.9)


