import sys
import os
sys.path.append('C:/Users/brake/Documents/pyNB')
sys.path.append('E:/Research_Projects/004_Propofol/Modelling/neuron_simulations')
import niktools as nt
import numpy as np
from matplotlib import pyplot as plt
import LFPy
import fig1
from neuron import h
import mat4py as m4p


np.random.seed(1)

###################
### Build model ###
###################

# file = 'E:/Research_Projects/Propofol/Data/dendrite_morphology/prinicple_cells/04a_pyramidal4aACC'
# file = 'E:/Research_Projects/Propofol/Data/dendrite_morphology/NMO_86952_layer 3_interneuron_Aspiny'
file = 'E:/Research_Projects/004_Propofol/Modelling/neuron_simulations/data/dendrite_morphologies/prinicple_cells/86955'
# file = 'E:/Research_Projects/Propofol/Data/dendrite_morphology/prinicple_cells/NMO_36031'
pts3d,connections,segs,morphData = fig1.morph2Segs(file)
dendIdcs = np.concatenate(segs)
p0 = np.mean(morphData[dendIdcs-1,2:5],axis=0)


l = list()
for i in range(len(segs)):
    s = np.array(morphData[segs[i]-1,3:6])
    l.append(np.sum(np.linalg.norm(np.diff(s,axis=0),2,1)))

L = np.sum(l)

Epas = -60 # a la Ahlfors and Wreh
Gpas = 1/5000
Ra = 100 

# Add soma
soma = h.Section(name='soma')
soma.insert('pas')
soma.e_pas = Epas
soma.g_pas = Gpas
soma.Ra = Ra

soma.push()
h.pt3dclear()
h.pt3dadd(morphData[1,2], morphData[1,3], morphData[1,4], morphData[1,5])
h.pt3dadd(morphData[0,2], morphData[0,3], morphData[0,4], morphData[0,5])
h.pt3dadd(morphData[2,2], morphData[2,3], morphData[2,4], morphData[2,5])
sectionList = h.SectionList()
sectionList.append()
h.pop_section()

soma.nseg = 5

# Add 3D points for dendritic tree
dend = list()
x = list()
for i in range(pts3d.shape[0]):
    dend.append(h.Section(name='dend'+str(i)))
    dend[-1].push()
    h.pt3dclear()
    h.pt3dadd(pts3d[i,0], pts3d[i,1], pts3d[i,2], pts3d[i,3])
    h.pt3dadd(pts3d[i,4], pts3d[i,5], pts3d[i,6], pts3d[i,7])
    dend[-1].insert('pas')
    dend[-1].e_pas = Epas
    # x.append(np.sqrt(0.25*(pts3d[i,0]+pts3d[i,4])**2+0.25*(pts3d[i,1]+pts3d[i,5])**2+0.25*(pts3d[i,2]+pts3d[i,6])**2)) # Distance to soma
    # dend[-1].g_pas = 1e-3/(8 + 44/(1+np.exp((x[-1]-406)/50))) # Destexhe & Pare (1999)
    dend[-1].g_pas = Gpas
    dend[-1].nseg = 3
    dend[-1].Ra = Ra
    sectionList.append()
    h.pop_section()

# Connect compartments together
for i in range(len(connections)):
    i1 = connections[i][0]
    i0 = connections[i][1]
    if(i1==-1):
        dend[i0].connect(soma,1,0)
    else:
        dend[i0].connect(dend[i1](1))


################
### Simulate ###
################

propofol = int(sys.argv[1])
tmax = 1000

cell = LFPy.Cell(sectionList,v_init=Epas,celsius=37)
cell.tstop = tmax
cell.simulate(rec_vmem=True)
t = cell.tvec

N = np.random.poisson(L*1.2)
# N = 100

N2 = len(t)
tmax = np.max(t)
N3 = cell.totnsegs

synIdcs = cell.get_rand_idx_area_norm(section='allsec', nidx=N, z_min=- 1000000.0, z_max=1000000.0)

def getQfromSynapse(synParams,time,sectionList):
    cell = LFPy.Cell(sectionList,v_init=Epas,celsius=37)
    cell.tstop = tmax
    syn = LFPy.Synapse(cell, **synParams)
    syn.set_spike_times(time)
    cdm = LFPy.CurrentDipoleMoment(cell=cell)
    cell.tstop = 100
    cell.simulate(probes=[cdm],rec_vmem= False)
    synPos0 = np.array([syn.x,syn.y,syn.z])
    return cdm, synPos0

synType = 'Exp2Syn'
iSyn = {'tau1':2,
        'tau2':7*(1+1.5*propofol),
        'weight':0.0001*(1+1.5*propofol),
        'erev':-80}
eSyn = {'tau1':0.5,
        'tau2':2,
        'weight':0.0002,
        'erev':0}

synParamsI = {'idx': 0,'e': iSyn['erev'],'syntype': synType,'tau1': iSyn['tau1'],'tau2': iSyn['tau2'],'weight': iSyn['weight'],'record_current': False}
synParamsE = {'idx': 0,'e': eSyn['erev'],'syntype': synType,'tau1': eSyn['tau1'],'tau2': eSyn['tau2'],'weight': eSyn['weight'],'record_current': False}

lambdaI = 1
lambdaE = 0.1



fE = lambda t,L: L+0.1*L*np.sin(2*np.pi*t*10/1000)
fI = lambda t,L: L+0.1*L*np.sin(2*np.pi*t*10/1000+np.pi)
tV = np.linspace(0,tmax,10000)
yI = np.cumsum(fI(tV,lambdaI))/np.sum(fI(tV,lambdaI))
yE = np.cumsum(fE(tV,lambdaE))/np.sum(fE(tV,lambdaE))

cell = LFPy.Cell(sectionList,v_init=Epas,celsius=37)
cell.tstop = tmax
syn = list()
synData = list()
for i in synIdcs:
    if(np.random.random()<=0.15):
        synParams = synParamsI.copy()
        n = np.random.poisson(lambdaI*tmax/1000)
        # ts = np.interp(np.random.random(n),yI,tV) # Periodic 10 Hz synaptic input
        ts = np.random.random(n)*tmax # Uniformly random synaptic input
    else:
        synParams = synParamsE.copy()
        n = np.random.poisson(lambdaE*tmax/1000)
        # ts = np.interp(np.random.random(n),yE,tV) # Periodic 10 Hz synaptic input
        ts = np.random.random(n)*tmax # Uniformly random synaptic input
    synParams['idx'] = i.tolist()
    syn.append(LFPy.Synapse(cell, **synParams))
    syn[-1].set_spike_times(ts)
    synParams['time'] = ts.tolist() 
    synParams['position'] = np.array([syn[-1].x,syn[-1].y,syn[-1].z]).tolist()
    synData.append(synParams)

cdm = LFPy.CurrentDipoleMoment(cell=cell)
cell.simulate(probes=[cdm],rec_vmem= True)
V = cell.vmem
Q = cdm.data
Qx = Q[0,:]
Qy = Q[1,:]
Qz = Q[2,:]

simulatenous = {'Qx':Qx.tolist(),
    'Qy':Qy.tolist(),
    'Qz':Qz.tolist(),
    'V':V.tolist()}

Q = np.zeros((N,3,N2))
synPos = np.zeros((N,3))
# for i in range(N):
#     print('Simulation ' + str(i+1) + ': ',end='')
#     synParams['idx'] = synIdcs[i]
#     cdm,synPos0 = getQfromSynapse(synParams,np.array(synTime[i]),sectionList)
#     Q[i,:,:] = cdm.data
#     synPos[i,:] = synPos0

Qx = Q[:,0,:]
Qy = Q[:,1,:]
Qz = Q[:,2,:]

independent = {'Qx':Qx.tolist(),
    'Qy':Qy.tolist(),
    'Qz':Qz.tolist()}

morphology = {'segments':[x.tolist() for x in segs],
            'coordinates':morphData.tolist(),
            'Vx':cell.x.tolist(),
            'Vy':cell.y.tolist(),
            'Vz':cell.z.tolist(),
            'fileName':file}

data = {'Epas':Epas,
    'Ra':Ra,
    't':t.tolist(),
    'synapses':synData,
    'morphology':morphology,
    'independent':independent,
    'simulatenous':simulatenous}



savePath = 'E:/Research_Projects/004_Propofol/Modelling/neuron_simulations/data/simulations/temp/'
nrn_name = os.path.split(file)[-1]
if(propofol):
    condition = '_propofol_simulation'
else:
    condition = '_simulation'


number = 1
savedName = savePath  + nrn_name + condition + f"{number:02d}" + '.mat'
while(os.path.exists(savedName)):
    number += 1
    savedName = savePath  + nrn_name + condition + f"{number:02d}" + '.mat'

m4p.savemat(savedName, data)

nt.bell()