from tensorflow.keras.utils import Sequence
import pandas as pd
import numpy as np
import random
import math
import pysam
from .util import *
import threading
import pickle
import pdb

def get_weights(data):
    w1=[float(data.shape[0])/sum(data.iloc[:,i]==1) for i in range(data.shape[1])]
    w0=[float(data.shape[0])/sum(data.iloc[:,i]==0) for i in range(data.shape[1])]
    return w1,w0


def get_probability_thresh_for_precision(truth,predictions,precision_thresh):
    from sklearn.metrics import precision_recall_curve
    num_tasks=truth.shape[1]
    precision_thresholds=[]
    for task_index in range(num_tasks):
        truth_task=truth.iloc[:,task_index]
        pred_task=predictions[:,task_index]
        non_ambig=truth_task!=-1
        precision,recall,threshold=precision_recall_curve(truth_task[non_ambig],pred_task[non_ambig])
        threshold=np.insert(threshold,threshold.shape[0],1)
        merged_prc=pd.DataFrame({'precision':precision,
                                 'recall':recall,
                                 'threshold':threshold})
        precision_thresholds.append(np.min(merged_prc[merged_prc['precision']>=precision_thresh]['threshold']))
    print(precision_thresholds)
    return precision_thresholds

def open_data_file(data_path=None,tasks=None,chroms_to_use=None):
    if data_path.endswith('.hdf5'):
        if tasks==None:
            data=pd.read_hdf(data_path)
        else:
            data=pd.read_hdf(data_path,columns=tasks)
    else:
        #treat as bed file
        if tasks==None:
            data=pd.read_csv(data_path,header=0,sep='\t')
        else:
            data=pd.read_csv(data_path,header=0,sep='\t',nrows=1)
            chrom_col=data.columns[0]
            start_col=data.columns[1]
            end_col=data.columns[2]
            data=pd.read_csv(data_path,header=0,sep='\t',usecols=[chrom_col,start_col,end_col]+tasks)
    print("loaded labels")
    try:
        data=data.set_index(['CHR','START','END'])
        print('set index to CHR, START, END')
    except:
        pass
    if chroms_to_use!=None:
        data=data[np.in1d(data.index.get_level_values(0), chroms_to_use)]
    print("filtered on chroms_to_use")
    print("data.shape:"+str(data.shape))
    return data

class TruePosGenerator(Sequence):
    def __init__(self,data_pickle,ref_fasta,batch_size=128,precision_thresh=0.9,expand_dims=True):
        f=open(data_pickle,'rb')
        data=pickle.load(f)
        self.predictions=data[0]
        self.labels=data[1]
        self.columns=self.labels.columns
        #calculate prediction probability cutoff to achieve the specified precision threshold
        self.prob_thresholds=get_probability_thresh_for_precision(self.labels,self.predictions,precision_thresh)
        truth_pred_product=self.labels*(self.predictions>=self.prob_thresholds)
        true_pos_rows=truth_pred_product[truth_pred_product.max(axis=1)>0]
        self.data=true_pos_rows
        self.indices=np.arange(self.data.shape[0])
        self.add_revcomp=False
        self.ref_fasta=ref_fasta
        self.lock=threading.Lock()
        self.batch_size=batch_size
        self.expand_dims=expand_dims


    def __len__(self):
        return math.ceil(self.data.shape[0]/self.batch_size)

    def __getitem__(self,idx):
        with self.lock:
            self.ref=pysam.FastaFile(self.ref_fasta)
            return self.get_basic_batch(idx)

    def get_basic_batch(self,idx):
        #get seq positions
        inds=self.indices[idx*self.batch_size:(idx+1)*self.batch_size]
        bed_entries=self.data.index[inds]
        #get sequences
        seqs=[self.ref.fetch(i[0],i[1],i[2]) for i in bed_entries]
        #one-hot-encode the fasta sequences
        seqs=np.array([[ltrdict.get(x,[0,0,0,0]) for x in seq] for seq in seqs])
        x_batch=seqs
        if (self.expand_dims==True):
            x_batch=np.expand_dims(x_batch,1)
        #extract the labels at the current batch of indices
        y_batch=np.asarray(self.data.iloc[inds])
        return (bed_entries,x_batch,y_batch)


#use wrappers for keras Sequence generator class to allow batch shuffling upon epoch end
class DataGenerator(Sequence):
    def __init__(self,
                 index_path,
                 input_path,
                 output_path,
                 num_inputs,
                 num_outputs,
                 ref_fasta=None,
                 batch_size=128,
                 add_revcomp=True,
                 index_tasks=None,
                 tasks=None,
                 shuffled_ref_negatives=False,
                 chroms_to_use=None,
                 get_w1_w0=False,
                 expand_dims=True,
                 upsample_thresh_list=None,
                 upsample_ratio_list=None,
                 shuffle=True):
        
        self.lock = threading.Lock()
        self.expand_dims=expand_dims
        self.shuffle=shuffle
        self.batch_size=batch_size
        self.ref_fasta=ref_fasta
        self.chroms_to_use=chroms_to_use
        
        #decide if reverse complement should be used
        self.add_revcomp=add_revcomp
        if add_revcomp==True:
            self.batch_size=int(batch_size/2)

        #determine whether negative set should consist of the shuffled refs.
        # If so, split batch size in 2, as each batch will be augmented with shuffled ref negatives
        # in ratio equal to positives
        self.shuffled_ref_negatives=shuffled_ref_negatives
        if self.shuffled_ref_negatives==True:
            self.batch_size=int(self.batch_size/2)

        #get the index, input, and output files 
        self.index_tasks=index_tasks
        self.tasks=tasks        
        self.index_path=index_path,
        self.input_path=input_path,
        self.output_path=output_path,
        self.num_inputs=num_inputs,
        self.num_outputs=num_outputs,
        self.file_to_pd=get_file_to_pd()        
        self.indices=self.file_to_pd[self.index_path]
        num_indices=self.indices.shape[0]

        #handle task-specific weights -- this is a bit outdated and may be removed in the future. 
        if get_w1_w0==True:
            assert self.data is not None
            w1,w0=get_weights(self.data)
            self.w1=w1
            self.w0=w0
            
        #set variables needed for upsampling the positives
        self.upsample_thresh_list=upsample_thresh_list
        self.upsample_ratio_list=upsample_ratio_list
        #generate the upsampled threshold index subgroups 
        if self.upsample_thresh_list is not None:
            self.get_upsampled_indices(num_indices)
        else:
            self.indices=self.indices.index
            if self.shuffle == True:
                np.random.shuffle(self.indices)

    def get_file_to_pd(self):
        '''
        make sure all input/output/index files are loaded only once in case there's overlap 
        generate a dictionary of file name to pandas data frame of file contents 
        '''
        file_to_df={}
        file_to_df[self.index_path]=open_data_file(data_path=self.index_path,tasks=self.index_tasks,chroms_to_use=self.chroms_to_use)
        for i in range(self.num_inputs):
            cur_input=self.input_path[i]
            if cur_input=="seq":
                continue
            if cur_input in file_to_df:
                continue
            file_to_df[self.input_path[i]]=open_data_file(data_path=cur_input,tasks=self.tasks[i],chroms_to_use=self.chroms_to_use)
        for i in range(self.num_outputs):
            cur_output=self.output_path[i]
            if cur_output in file_to_df:
                continue
            file_to_df[cur_output]=open_data_file(data_path=cur_output,tasks=self.tasks[i],chroms_to_use=self.chroms_to_use)
        return file_to_df
            
    def get_upsampled_indices(self,num_indices):
        '''
        several levels of upsampling are handled
        self.upsample_thresh_list -- list of thresholds for upsampling the dataset 
        self.upsample_ratio_list -- fraction of batch to be generated at each threshold 
        '''

        self.upsampled_coord_indices = {}
        self.upsampled_numerical_indices = {}
        self.batch_sizes = []
        self.wraps = []

        print("upsample thresh list:"+str(self.upsample_thresh_list))
        if self.upsample_ratio_list is None:
            #all thresholds represented equally in the batch 
            self.upsample_ratio = 1 / (len(self.upsample_thresh_list) - 1)
            self.upsample_ratio_list = [self.upsample_ratio for i in range(len(self.upsample_thresh_list) - 1)]
            print("upsample ratio list: " + str(self.upsample_ratio_list))
            
        for ind,val in enumerate(self.upsample_thresh_list):
            #iterate through upsampling thresholds 
            if ind < (len(self.upsample_thresh_list) - 1):
                if ind == (len(self.upsample_thresh_list) - 2):
                    #determine the  number of examples from the batch size that should be selected from this threshold group 
                    self.batch_sizes.append(int(self.batch_size - sum(self.batch_sizes)))
                else:
                    self.batch_sizes.append(int(self.batch_size * self.upsample_ratio_list[ind]))
                #identify the coordinates that fall into the current threshold group
                self.upsampled_coord_indices[ind] = self.indices.loc[(self.indices >= val).any(axis=1) & (self.indices < self.upsample_thresh_list[ind+1]).any(axis=1)].index
                self.upsampled_numerical_indices[ind] = np.arange(self.upsampled_coord_indices[ind].shape[0])
                self.wraps.append(math.ceil(num_indices/self.upsampled_numerical_indices[ind].shape[0]))
                self.upsampled_numerical_indices[ind] = np.repeat(self.upsampled_numerical_indices[ind], self.wraps[ind])[0:num_indices]
                if self.shuffle == True:
                    np.random.shuffle(self.upsampled_numerical_indices[ind])
        return
            
        
    def __len__(self):
        return math.ceil(len(self.indices)/self.batch_size)

    def get_coords(self,idx):
        if self.upsample_thresh_list is not None:
            all_bed_entries=[]
            for ind,val in enumerate(self.batch_sizes):
                batch_indices = self.upsample_indices[ind][idx*val:(idx+1)*val]
                bed_entries = self.upsample_data[ind].index[batch_indices]
                all_bed_entries+=bed_entries
        else:
            inds=self.indices[idx*self.batch_size:(idx+1)*self.batch_size]
            all_bed_entries=self.indices[inds]
        return all_bed_entries
    
    def get_seq(self,coords):
        pass
    
    def get_pd_vals(self,coords):
        pass
    
    def transform_seq(self,seq):
        pass
    
    def one_hot_encode(self,seq):
        pass
    
    def transform_vals(self,vals):
        pass
    
    def __getitem__(self,idx):
        with self.lock:
            ref=pysam.FastaFile(self.ref_fasta)
            self.ref=ref
            
            #get the coordinates for the current batch
            coords=self.get_coords(idx)

            #get the inputs
            X=[]
            for cur_input_index in range(self.num_inputs):
                cur_input=self.input_path[cur_input_index]
                if cur_input=="seq":
                    X.append(self.one_hot_encode(self.transform_seq(self.get_seq(coords))))
                else:
                    #extract values from pandas df 
                    X.append(self.transform_vals(self.get_pd_vals(coords)))
            #get the outputs
            y=[]
            for cur_output_index in range(self.num_outputs):
                cur_output=self.output_path[cur_output_index] 
                if cur_ouput=="seq":
                    y.append(self.one_hot_encode(self.transform_seq(self.get_seq(coords)))
                else:
                    y.append(self.tranform_vals(self.get_pd_vals(coords)))                    
            #return the batch as an X,y tuple 
            return (X,y)
        
            '''
            if self.shuffled_ref_negatives==True:
                (X_seq,y,bed_entries)=self.get_shuffled_ref_negatives_batch(idx)
            elif self.upsample_ratio > 0:
                (X_seq,y,bed_entries)=self.get_upsampled_positives_batch(idx)
            elif self.upsample_ratio_list > 0:
                (X_seq,y,bed_entries)=self.get_upsampled_positives_batch(idx)
            else:
                (X_seq,y,bed_entries)=self.get_basic_batch(idx)
            ''' 
                          
    def get_shuffled_ref_negatives_batch(self,idx):
        #get seq positions
        inds=self.indices[idx*self.batch_size:(idx+1)*self.batch_size]
        bed_entries=self.data.index[inds]
        #get sequences
        seqs=[self.ref.fetch(i[0],i[1],i[2]) for i in bed_entries]
        if self.add_revcomp==True:
            #add in the reverse-complemented sequences for training.
            seqs_rc=[revcomp(s) for s in seqs]
            seqs=seqs+seqs_rc
            bed_entries=bed_entries+bed_entries
            
        #generate the corresponding negative set by dinucleotide-shuffling the sequences
        seqs_shuffled=[dinuc_shuffle(s) for s in seqs]
        seqs=seqs+seqs_shuffled
        bed_entries=bed_entries+bed_entries
        
        #one-hot-encode the fasta sequences
        seqs=np.array([[ltrdict.get(x,[0,0,0,0]) for x in seq] for seq in seqs])
        x_batch=seqs
        if (self.expand_dims==True):
            x_batch=np.expand_dims(x_batch,1)
        y_batch=np.asarray(self.data.iloc[inds])
        if self.add_revcomp==True:
            y_batch=np.concatenate((y_batch,y_batch),axis=0)
        y_shape=y_batch.shape
        y_batch=np.concatenate((y_batch,np.zeros(y_shape)))
        return (x_batch,y_batch,bed_entries)

    def get_upsampled_positives_batch(self,idx):
        #get seq positions
        if self.upsample_thresh_list is not None:
            seqs = []
            labels = []
            all_bed_entries=[]
            for ind,val in enumerate(self.batch_sizes):
                batch_indices = self.upsample_indices[ind][idx*val:(idx+1)*val]
                bed_entries = self.upsample_data[ind].index[batch_indices]
                seq_entries = [self.ref.fetch(i[0],i[1],i[2]) for i in bed_entries]
                seqs += seq_entries
                all_bed_entries+=bed_entries
                if ind == 0:
                    labels = self.upsample_data[ind].iloc[batch_indices].values
                else:
                    labels = np.append(labels,self.upsample_data[ind].iloc[batch_indices].values)
            y_batch = np.array(labels)
            y_batch = np.expand_dims(y_batch,1)

        if self.upsample_ratio > 0:
            pos_inds=self.pos_indices[idx*self.pos_batch_size:(idx+1)*self.pos_batch_size]
            pos_bed_entries=self.ones.index[pos_inds]
            neg_inds=self.neg_indices[idx*self.neg_batch_size:(idx+1)*self.neg_batch_size]
            try:
                neg_bed_entries=self.zeros.index[neg_inds]
            except:
                neg_bed_entries=self.zeros[neg_inds]

            #get sequences
            pos_seqs=[self.ref.fetch(i[0],i[1],i[2]) for i in pos_bed_entries]
            neg_seqs=[self.ref.fetch(i[0],i[1],i[2]) for i in neg_bed_entries]
            seqs=neg_seqs+pos_seqs

            bed_entries=[(i[0],i[1],i[2]) for i in pos_bed_entries]+[(i[0],i[1],i[2]) for i in neg_bed_entries]
            
        if self.add_revcomp==True:
            #add in the reverse-complemented sequences for training.
            seqs_rc=[revcomp(s) for s in seqs]
            seqs=seqs+seqs_rc
            bed_entries+=bed_entries

        #one-hot-encode the fasta sequences
        seqs=np.array([[ltrdict.get(x,[0,0,0,0]) for x in seq] for seq in seqs])
        x_batch=seqs

        if (self.expand_dims==True):
            x_batch=np.expand_dims(x_batch,1)

        if self.upsample_ratio > 0:
            #extract the positive and negative labels at the current batch of indices
            y_batch_pos=self.ones.iloc[pos_inds]
            y_batch_neg=self.zeros.iloc[neg_inds]
            y_batch=np.concatenate((y_batch_neg,y_batch_pos),axis=0)
            
        #add in the labels for the reverse complement sequences, if used
        if self.add_revcomp==True:
            y_batch=np.concatenate((y_batch,y_batch),axis=0)
        return (x_batch,y_batch,bed_entries)

    def get_basic_batch(self,idx):
        #get seq positions
        inds=self.indices[idx*self.batch_size:(idx+1)*self.batch_size]
        if self.data_path is not None:
            bed_entries=self.data.index[inds]
        else:
            pos_inds=inds[inds<self.universal_negative_offset]
            neg_inds=inds[inds>=self.universal_negative_offset]
            bed_entries=np.concatenate((self.nonzero_bins.index[pos_inds],self.universal_negatives.index[neg_inds-self.universal_negative_offset]),axis=0)
        #get sequences
        seqs=[self.ref.fetch(i[0],i[1],i[2]) for i in bed_entries]
        if self.add_revcomp==True:
            #add in the reverse-complemented sequences for training.
            seqs_rc=[revcomp(s) for s in seqs]
            seqs=seqs+seqs_rc
            bed_entries+=bed_entries
            
        #one-hot-encode the fasta sequences
        seqs=np.array([[ltrdict.get(x,[0,0,0,0]) for x in seq] for seq in seqs])
        x_batch=seqs
        if(self.expand_dims==True):
            x_batch=np.expand_dims(x_batch,1)
        #extract the labels at the current batch of indices
        if self.data_path is not None:
            y_batch=np.asarray(self.data.iloc[inds])
        else:
            y_batch_pos=self.nonzero_bins.iloc[pos_inds]
            y_batch_neg=np.zeros((x_batch.shape[0]-y_batch_pos.shape[0],y_batch_pos.shape[1]))
            y_batch=np.concatenate((y_batch_pos,y_batch_neg),axis=0)
        #add in the labels for the reverse complement sequences, if used
        if self.add_revcomp==True:
            y_batch=np.concatenate((y_batch,y_batch),axis=0)
        return (x_batch,y_batch,bed_entries)

    def on_epoch_end(self):
        #if upsampling is being used, shuffle the positive and negative indices
        if self.shuffle==True:
            if self.upsample_thresh_list is not None:
                for ind,val in enumerate(self.batch_sizes):
                    np.random.shuffle(self.upsample_indices[ind])
            elif self.upsample_ratio > 0:
                np.random.shuffle(self.pos_indices)
                np.random.shuffle(self.neg_indices)
            else:
                np.random.shuffle(self.indices)


#generate batches of SNP data with specified allele column name and flank size
class SNPGenerator(DataGenerator):
    def __init__(self,
                 allele_col,
                 flank_size,
                 data_path=None,
                 nonzero_bin_path=None,
                 universal_negative_path=None,
                 ref_fasta=None,
                 batch_size=128,
                 add_revcomp=True,
                 tasks=None,
                 shuffled_ref_negatives=False,
                 upsample_thresh=0,
                 upsample_ratio=0,
                 chroms_to_use=None,
                 get_w1_w0=False,
                 expand_dims=True,
                 upsample_thresh_list=None,
                 upsample_ratio_list=None,
                 shuffle=True):

        DataGenerator.__init__(self,
                               data_path=data_path,
                               nonzero_bin_path=nonzero_bin_path,
                               universal_negative_path=universal_negative_path,
                               ref_fasta=ref_fasta,
                               batch_size=batch_size,
                               add_revcomp=add_revcomp,
                               tasks=tasks,
                               shuffled_ref_negatives=shuffled_ref_negatives,
                               upsample_thresh=upsample_thresh,
                               upsample_ratio=upsample_ratio,
                               chroms_to_use=chroms_to_use,
                               get_w1_w0=get_w1_w0,
                               expand_dims=expand_dims,
                               upsample_thresh_list=upsample_thresh_list,
                               upsample_ratio_list=upsample_ratio_list,
                               shuffle=shuffle)
        
        self.allele_col=allele_col
        self.flank_size=flank_size

    #override the get_basic_batch definition to extract flanks around the bQTL
    # and insert the specified allele
    def get_basic_batch(self,idx):
        #get seq positions
        inds=self.indices[idx*self.batch_size:(idx+1)*self.batch_size]
        entries=self.data.iloc[inds]
        seqs=[]
        for index,row in entries.iterrows():
            allele=row[self.allele_col].split(',')[0]
            chrom=index[0]
            pos=index[1]-1
            start=pos-self.flank_size
            end=pos+self.flank_size
            seq=self.ref.fetch(chrom,start,end)
            left=seq[0:self.flank_size]
            right=seq[self.flank_size+len(allele)::]
            seq=seq[0:self.flank_size]+allele+seq[self.flank_size+len(allele)::]
            seqs.append(seq)
        #one-hot-encode the fasta sequences
        seqs=np.array([[ltrdict.get(x,[0,0,0,0]) for x in seq] for seq in seqs])
        x_batch=seqs
        if (self.expand_dims==True):
            x_batch=np.expand_dims(x_batch,1)
        return x_batch
