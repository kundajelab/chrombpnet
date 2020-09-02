rm(list=ls())
library(ggplot2)
source("~/helpers.R")
jsd=read.table("k562.dnase.jsd.txt",header=TRUE,sep='\t')
mnl=read.table("k562.dnase.mnnll.txt",header=TRUE,sep='\t')
p1=ggplot(data=jsd)+
  stat_ecdf(aes(x=jsd$baseline,color='baseline'))+
  stat_ecdf(aes(x=jsd$bpnet,color='BPNET'))+
  stat_ecdf(aes(x=jsd$tobias,color="Tobias"))+
  stat_ecdf(aes(x=jsd$bpnet.tobias,color='BPNET.Tobias'))+
  scale_color_manual(values=c("baseline"="#FF0000",
                              "BPNET"="#000000",
                              "Tobias"="#0000FF",
                              "BPNET.Tobias"="#00FF00"))+
  xlab("JSD Test Set")+
  ylab("CDF")+
  ggtitle("K562 DNASE")
p2=ggplot(data=mnl)+
  stat_ecdf(aes(x=mnl$baseline,color='baseline'))+
  stat_ecdf(aes(x=mnl$bpnet,color='BPNET'))+
  stat_ecdf(aes(x=mnl$tobias,color="Tobias"))+
  stat_ecdf(aes(x=mnl$bpnet.tobias,color='BPNET.Tobias'))+
  scale_color_manual(values=c("baseline"="#FF0000",
                              "BPNET"="#000000",
                              "Tobias"="#0000FF",
                              "BPNET.Tobias"="#00FF00"))+
  xlab("Multinomial log likelihood Test Set")+
  ylab("CDF")
multiplot(p1,p2,cols=1)

jsd=read.table("hepg2.dnase.jsd.txt",header=TRUE,sep='\t')
mnl=read.table("hepg2.dnase.mnnll.txt",header=TRUE,sep='\t')
p3=ggplot(data=jsd)+
  stat_ecdf(aes(x=jsd$baseline,color='baseline'))+
  stat_ecdf(aes(x=jsd$bpnet,color='BPNET'))+
  stat_ecdf(aes(x=jsd$tobias,color="Tobias"))+
  stat_ecdf(aes(x=jsd$bpnet.tobias,color='BPNET.Tobias'))+
  scale_color_manual(values=c("baseline"="#FF0000",
                              "BPNET"="#000000",
                              "Tobias"="#0000FF",
                              "BPNET.Tobias"="#00FF00"))+
  xlab("JSD Test Set")+
  ylab("CDF")+
  ggtitle("HEPG2 DNASE")
p4=ggplot(data=mnl)+
  stat_ecdf(aes(x=mnl$baseline,color='baseline'))+
  stat_ecdf(aes(x=mnl$bpnet,color='BPNET'))+
  stat_ecdf(aes(x=mnl$tobias,color="Tobias"))+
  stat_ecdf(aes(x=mnl$bpnet.tobias,color='BPNET.Tobias'))+
  scale_color_manual(values=c("baseline"="#FF0000",
                              "BPNET"="#000000",
                              "Tobias"="#0000FF",
                              "BPNET.Tobias"="#00FF00"))+
  xlab("Multinomial log likelihood Test Set")+
  ylab("CDF")
multiplot(p3,p4,cols=1)