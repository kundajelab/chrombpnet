import h5py
import numpy as np

# Load function
def load_dict_from_h5(h5group):
    dic = {}
    for key, item in h5group.items():
        if isinstance(item, h5py.Group):
            # If subgroup looks like a list (all keys are integers)
            if all(k.isdigit() for k in item.keys()):
                lst = [item[str(i)][()] for i in sorted(map(int, item.keys()))]
                dic[key] = lst
            else:
                dic[key] = load_dict_from_h5(item)
        elif isinstance(item, h5py.Dataset):
            dic[key] = item[()]
    return dic

# Function to load the whole dictionary from HDF5 file
def load_from_h5(filename):
    with h5py.File(filename, 'r') as h5file:
        return load_dict_from_h5(h5file)

loaded_data = load_from_h5('composite_footprint_data.h5')

print(loaded_data['HEPG2_profile-neg.pattern_4'].keys())
print("num distances",len(loaded_data['HEPG2_profile-neg.pattern_4']['combine_counts'].keys()))

print(loaded_data['HEPG2_profile-neg.pattern_4']['combine_counts'].keys())
print(loaded_data['HEPG2_profile-neg.pattern_4']['combine_footprint'].keys())

print(loaded_data['HEPG2_profile-neg.pattern_4']['control'].keys())
print(loaded_data['HEPG2_profile-neg.pattern_4']['motif1'].keys())
print(loaded_data['HEPG2_profile-neg.pattern_4']['motif2'].keys())

print(loaded_data['HEPG2_profile-neg.pattern_4']['motif1']['counts'])
print(loaded_data['HEPG2_profile-neg.pattern_4']['motif1']['profile'].shape)

print(loaded_data['HEPG2_profile-neg.pattern_4']['combine_counts']["5"])
print(loaded_data['GM12878_profile-neg.pattern_35']['combine_counts']["5"])




