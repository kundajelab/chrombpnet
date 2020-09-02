library(ggplot2)
source("~/helpers.R")
jsd=read.table("hepg2.dnase.jsd.txt",header=TRUE,sep='\t')
mnl=read.table("hepg2.dnase.mnnll.txt",header=TRUE,sep='\t')
jsd$delta=jsd$bpnet.tobias-jsd$baseline

p1=ggplot(data=jsd,aes(x=jsd$delta))+
  geom_histogram(bins=100)+
  xlab("HEPG2 DNASE Corrected BPNET.Tobias - Baseline")+
  ylab("N Test Set Examples")+
  ggtitle("JSD")

mnl$delta=mnl$bpnet.tobias-mnl$baseline
p2=ggplot(data=mnl,aes(x=mnl$delta))+
  geom_histogram(bins=100)+
  xlab("HEPG2 DNASE Corrected BPNET.Tobias - Baseline")+
  ylab("N Test Set Examples")+
  ggtitle("MNLL")
multiplot(p1,p2,cols=1)
merged=merge(jsd,mnl,by=c("chrom","summit"))
merged=merged[order(merged$delta.x),]
sub=subset(merged,select=c("chrom","summit","delta.x","delta.y"))
colnames(sub)=c("Chrom","Summit","deltaJSD","deltaMNLL")