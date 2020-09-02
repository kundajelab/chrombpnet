library(ggplot2)
args = commandArgs(trailingOnly=TRUE)
seed=args[1]
title=args[2] 
prefix=args[3] 
outf=args[4]
assay=args[6]

counts_loss_scale_factor=as.numeric(args[5])
data_batch=read.table(paste(prefix,"/seed",seed,"/bias.",assay,".0.log.batch",sep=''),header=TRUE,sep='\t')
data_epoch=read.table(paste(prefix,"/seed",seed,"/bias.",assay,".0.log",sep=''),header=TRUE,sep=',')
data_batch$logcount_predictions_loss_scaled=data_batch$logcount_predictions_loss*counts_loss_scale_factor
data_batch$Batch=seq(1,nrow(data_batch))

p1=ggplot(data=data_batch)+
  geom_line(aes(x=Batch,y=loss,color="Training_Total"))+
  geom_line(aes(x=Batch,y=logcount_predictions_loss,color="Training_LogCount"))+
  geom_line(aes(x=Batch,y=logcount_predictions_loss_scaled,color=paste("Training_LogCount_x",counts_loss_scale_factor,sep='')))+
  geom_line(aes(x=Batch,y=profile_predictions_loss,color="Training_Profile"))+
  xlab("Batch")+
  ylab("Loss")+
  ggtitle(paste(title,"Batch-level Loss"))+
  scale_color_manual(name="Loss",
                     values = c(
                       Training_LogCount='#1b9e77',
                       Training_LogCount_x18='#000000',
                       Training_Profile='#d95f02',
                       Training_Total='#7570b3'))+  
  scale_y_continuous(trans = 'log2')+
  theme_bw(20)

pdf(file=paste("batch",outf,sep=''),width=10,height=8,pointsize=12)
p1
dev.off()

data_epoch$logcount_predictions_loss_scaled=data_epoch$logcount_predictions_loss*counts_loss_scale_factor
data_epoch$val_logcount_predictions_loss_scaled=data_epoch$val_logcount_predictions_loss*counts_loss_scale_factor

p2=ggplot(data=data_epoch)+
  geom_line(aes(x=epoch,y=loss,color="Training_Loss"))+
  geom_line(aes(x=epoch,y=logcount_predictions_loss,color="Training_Counts_Loss"))+
  geom_line(aes(x=epoch,y=logcount_predictions_loss_scaled,color=paste("Training_Counts_Loss_x",counts_loss_scale_factor,sep="")))+
  geom_line(aes(x=epoch,y=profile_predictions_loss,color="Training_Profile_Loss"))+
  geom_line(aes(x=epoch,y=val_loss,color="Validation_Loss"))+
  geom_line(aes(x=epoch,y=val_logcount_predictions_loss,color="Validation_Counts_Loss"))+
  geom_line(aes(x=epoch,y=val_logcount_predictions_loss_scaled,color=paste("Validation_Counts_Loss_x",counts_loss_scale_factor,sep="")))+
  geom_line(aes(x=epoch,y=val_profile_predictions_loss,color="Validation_Profile_Loss"))+
  scale_color_manual(name="Loss",
                     values=c(Training_Loss="#e41a1c",
                              Training_Counts_Loss="#377eb8",
                              Training_Counts_Loss_x18="#a65628",
                              Training_Profile_Loss="#4daf4a",
                              Validation_Loss="#984ea3",
                              Validation_Counts_Loss_x18="#f781bf",
                              Validation_Counts_Loss="#ff7f00",
                              Validation_Profile_Loss="#ffff33"))+
  scale_y_continuous(trans='log2')+
  theme_bw(20)+
  ggtitle(paste(title,"Epoch Loss"))

pdf(file=paste("epoch",outf,sep=''),width=10,height=8,pointsize=12)
p2
dev.off()

