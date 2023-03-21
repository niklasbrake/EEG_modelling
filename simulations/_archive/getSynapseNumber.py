import sys
import os
sys.path.append('C:/Users/brake/Documents/pyNB')
sys.path.append('E:/Research_Projects/004_Propofol/Modelling/neuron_simulations')
sys.path.append('E:/Research_Projects/004_Propofol/Modelling/neuron_simulations/code')
import niktools as nt
import numpy as np
from matplotlib import pyplot as plt
import LFPy
import getMorphoSegments
from neuron import h, load_mechanisms
import mat4py as m4p
from scipy.stats import poisson


filepath = 'E:/Research_Projects/004_Propofol/Modelling/neuron_simulations/data/cortical_column_Hagen/swc/'
neurons = os.listdir(filepath)

L = list()
for neuron in neurons:
    file = filepath + neuron
    pts3d,connections,segs,morphData = getMorphoSegments.morph2Segs(file)
    dendIdcs = np.concatenate(segs)
    p0 = np.mean(morphData[dendIdcs-1,2:5],axis=0)
    l = list()
    for i in range(len(segs)):
        s = np.array(morphData[segs[i]-1,3:6])
        l.append(np.sum(np.linalg.norm(np.diff(s,axis=0),2,1)))
    L.append(1.15*np.sum(l))

print(L)