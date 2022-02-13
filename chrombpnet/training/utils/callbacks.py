import tensorflow.keras as keras

class LossHistory(keras.callbacks.Callback):
    """
    Callbacks to store train, validation loss at the the end of every batch and the end of every epoch.
    You can also track the counts loss and profile loss seperatley using the callbacks provided.
    """
    
    def __init__(self,model_output_path_logs_name,to_track):
        self.model_output_path_logs_name=model_output_path_logs_name
        self.to_track=to_track
        self.outf=open(self.model_output_path_logs_name,'w')
        self.outf.write('Epoch\tBatch\t'+'\t'.join(self.to_track)+'\n')
        keras.callbacks.Callback.__init__(self)

    def on_train_begin(self, logs={}):
        self.losses ={}

    def on_epoch_begin(self,epoch, logs={}):
        self.losses[epoch]={}
        for trackable in self.to_track:
            self.losses[epoch][trackable]=[]         
        self.cur_epoch=epoch
        
    def on_batch_end(self, batch, logs={}):
        for trackable in self.to_track:
            self.losses[self.cur_epoch][trackable].append(logs.get(trackable))
        
    def on_epoch_end(self,epoch,logs={}):
        marker=self.to_track[0] 
        num_batches=len(self.losses[self.cur_epoch][marker])
        for i in range(num_batches):
            self.outf.write(str(epoch)+'\t'+str(i))
            for trackable in self.to_track:
                self.outf.write('\t'+str(self.losses[epoch][trackable][i]))
            self.outf.write('\n')

        
    def on_train_end(self,logs={}):
        self.outf.close()
