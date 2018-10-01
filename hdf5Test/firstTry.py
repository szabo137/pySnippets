"""
first try of h5py usage - writer
"""
import h5py
import numpy as np
import time
import os

arrR = np.random.randn(1000)
arrI = 1j*np.random.randn(1000)
arr = arrR + arrI

with h5py.File('random.hdf5','w') as f:
    dset = f.create_dataset("default",data=arr)


baseMeta = {
    'Date': time.time(),
    'User':'me',
    'OS':os.name
}

with h5py.File('groups.hdf5','w') as f:
    g = f.create_group("base_complex")
    g.attrs.update(baseMeta)
    
    gr=g.create_group("sub_real")
    gi=g.create_group("sub_imag")
    
    d = g.create_dataset('default',data=arr)
    dr=gr.create_dataset('default',data=arrR)
    di=gi.create_dataset('default',data=arrI)


print "json -----------"
import json
jmeta = json.dumps(baseMeta)
print jmeta
print type(jmeta)

recoverMeta = json.loads(jmeta)
print recoverMeta
print type(recoverMeta)


print "pickle ------------"
import cPickle as pickle
pmeta = pickle.dumps(baseMeta)

print pmeta
print type(pmeta)

recoverMeta2 = pickle.loads(pmeta)
print recoverMeta2
print type(recoverMeta2)
