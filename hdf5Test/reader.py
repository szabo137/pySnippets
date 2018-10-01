"""
first try of h5py usage - reader 
"""
import h5py
import numpy as np


print "reading data ------------------------"
with h5py.File('random.hdf5','r') as f:
    data = f['default']
    print type(data)
    #print dir(data)
    print max(data)
    print data[:15]


print "selective reading --------------------"
with h5py.File('random.hdf5','r') as f:
    data_set = f['default']
    data = data_set[:10]
print data

print "group reading -------------------------"
with h5py.File("groups.hdf5",'r') as f:
    d=f['base_complex/default']
    dr=f["base_complex/sub_real/default"]
    di=f["base_complex/sub_imag/default"]
    print d[:10]
    print dr[:10]
    print di[:10]
    
print "group reading (nested loop)-------------------------"
with h5py.File("groups.hdf5",'r') as f:
    for k in f.keys():
        print "-"*20
        print k
        for kk in f[k].keys():
            print "\t%s"%kk
        

print "group reading (visit)-------------------------"

def get_all(name):
    if 'imag' in name:
        return name

with h5py.File("groups.hdf5",'r') as f:
    g= f.visit(get_all)
    print g
    

print "group reading (attrs)-------------------------"

with h5py.File("groups.hdf5",'r') as f:
    g=f['base_complex']
    for k in g.attrs.keys():
        print "%s: %s"%(k,g.attrs[k])
