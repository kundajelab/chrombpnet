import math 
import tiledb

def transform_data_type_min(inputs,num_inputs):
    if inputs is None:
        transformed=[-math.inf]*num_inputs
    elif inputs is []:
        transformed=[-math.inf]*num_inputs
    else:
        assert(len(inputs)==num_inputs)
        transformed=[]
        for i in range(num_inputs):
            transformed.append([]) 
            cur_inputs=inputs[i].split(',')
            for j in cur_inputs: 
                if str(j).lower()=="none":
                    transformed[i].append(-math.inf)
                else:
                    transformed[i].append(float(j))
    return transformed

def transform_data_type_max(inputs,num_inputs):
    if inputs is None:
        transformed=[math.inf]*num_inputs
    elif inputs is []:
        transformed=[math.inf]*num_inputs
    else:
        assert(len(inputs)==num_inputs)
        transformed=[]
        for i in range(num_inputs):
            transformed.append([]) 
            cur_inputs=inputs[i].split(',')
            for j in cur_inputs: 
                if str(j).lower()=="none":
                    transformed[i].append(math.inf)
                else:
                    transformed[i].append(float(j))
    return transformed

## not sure how important these are check with Anna - ask Anna to clean up

tdb_config_params={"sm.check_coord_dups":False,
                   "sm.check_coord_oob":False,
                   "sm.check_global_order":False,
                   "sm.num_writer_threads":40,
                   "sm.num_reader_threads":40,
                   "sm.num_async_threads":40,
                   "vfs.num_threads":40}

def get_default_config():
    tdb_config=tiledb.Config()
    tdb_config['vfs.s3.region']='us-west-1'
    tdb_config["sm.check_coord_dups"]="false"
    tdb_config["sm.check_coord_oob"]="false"
    tdb_config["sm.check_global_order"]="false"
    tdb_config["sm.num_reader_threads"]="40"
    tdb_config["sm.num_async_threads"]="40"
    tdb_config["vfs.num_threads"]="40"
    tdb_config["vfs.file.enable_filelocks"]="false"
    return tdb_config

