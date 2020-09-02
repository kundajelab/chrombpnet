rm(list=ls())
library(ggplot2)
source("~/helpers.R")
jsd=read.table("gm12878.atac.jsd.txt",header=TRUE,sep='\t')
mnl=read.table("gm12878.atac.mnnll.txt",header=TRUE,sep='\t')
p1=ggplot(data=jsd)+
  stat_ecdf(aes(x=jsd$baseline,color='baseline'))+
  stat_ecdf(aes(x=jsd$tobias,color="Tobias"))+
  scale_color_manual(values=c("baseline"="#FF0000",
                              "Tobias"="#0000FF"))+
  xlab("JSD Test Set")+
  ylab("CDF")+
  ggtitle("GM12878 ATAC")
p2=ggplot(data=mnl)+
  stat_ecdf(aes(x=mnl$baseline,color='baseline'))+
  stat_ecdf(aes(x=mnl$tobias,color="Tobias"))+
  scale_color_manual(values=c("baseline"="#FF0000",
                              "Tobias"="#0000FF"))+
  xlab("Multinomial log likelihood Test Set")+
  ylab("CDF")
multiplot(p1,p2,cols=1)


jsd=read.table("hepg2.atac.jsd.txt",header=TRUE,sep='\t')
mnl=read.table("hepg2.atac.mnnll.txt",header=TRUE,sep='\t')
p3=ggplot(data=jsd)+
  stat_ecdf(aes(x=jsd$baseline,color='baseline'))+
  stat_ecdf(aes(x=jsd$tobias,color="Tobias"))+
  scale_color_manual(values=c("baseline"="#FF0000",
                              "Tobias"="#0000FF"))+
  xlab("JSD Test Set")+
  ylab("CDF")+
  ggtitle("HEPG2 ATAC")
p4=ggplot(data=mnl)+
  stat_ecdf(aes(x=mnl$baseline,color='baseline'))+
  stat_ecdf(aes(x=mnl$tobias,color="Tobias"))+
  scale_color_manual(values=c("baseline"="#FF0000",
                              "Tobias"="#0000FF"))+
  xlab("Multinomial log likelihood Test Set")+
  ylab("CDF")
multiplot(p3,p4,cols=1)