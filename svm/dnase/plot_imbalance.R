rm(list=ls())
library(ggplot2)
library(dplyr)
data=read.table("genomewide.imbalance.txt",header=TRUE,sep='\t')
data$expectedAUPRC=data$NumPositives/(data$NumPositives+data$NumNegatives)
data$Fold=factor(data$Fold,levels=seq(0,9))
fold_means=data %>%
  group_by(Fold) %>%
  dplyr::summarize(Mean = mean(expectedAUPRC, na.rm=TRUE))
fold_means=fold_means[order(fold_means$Mean),]
data$Fold=factor(data$Fold,levels=fold_means$Fold)
attach(data)
p1=ggplot(data, aes(x = Fold,y=expectedAUPRC,fill=Task)) +
  geom_bar(stat='identity',position='dodge') +
  scale_fill_manual(values=c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00'))+
  theme_bw(15)+
  ggtitle("Canonical Cell Lines Genomewide expected auPRC")


