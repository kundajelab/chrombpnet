rm(list=ls())
library(ggplot2)
data=read.table("aggregate_performance_for_plot.tsv",header=TRUE,sep='\t')
# columns to paste together
cols <- c( 'Model' , 'TrainSet' , 'TestSet' )

# create a new column `x` with the three columns collapsed together
data$Experiment <- apply( data[ , cols ] , 1 , paste , collapse = "-" )
data$Experiment=factor(data$Experiment,levels=c("SVM-SVM-Genome",
                                                "CNN-SVM-Genome",
                                                "CNN-Genome-Genome",
                                                "SVM-SVM-SVM",
                                                "CNN-Genome-SVM",
                                                "CNN-SVM-SVM"))
data=data[order(data$Experiment),]
attach(data)
p1=ggplot(data, aes(x = Task,y=auPRC,fill=Task)) +
  geom_boxplot() +
  geom_jitter(color='black',shape = 21)+
  facet_grid(.~Experiment)+
  scale_fill_manual(values=c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00'))+
  theme_bw(15)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ggtitle("Model - Training Set - Test Set")
  

  