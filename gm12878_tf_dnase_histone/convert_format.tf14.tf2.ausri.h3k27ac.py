import h5py
#f = h5py.File("gm12878.dnase.with.bias.unplugged.0.hdf5",'r+')
#f = h5py.File("Spi1ProfileModel_RevComp_9layers.h5",'r+')
#f = h5py.File("patience_6_5K_out_9_gm12878.h3k27ac.seed.2345.cs.25.filters.300.naive.range.4.6.to.11.5.0.hdf5",'r+')
f = h5py.File("patience_6_5K_out_12_gm12878.h3k27ac.seed.2345.cs.25.filters.300.naive.range.4.6.to.11.5.0.hdf5",'r+')
data_p = f.attrs['training_config']
data_p = data_p.decode().replace("learning_rate","lr").encode()
f.attrs['training_config'] = data_p
f.close()
