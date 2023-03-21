import sys
import os
sys.path.append('C:/Users/brake/Documents/pyNB')
sys.path.append('E:/Research_Projects/004_Propofol/Modelling/neuron_simulations')
sys.path.append('E:/Research_Projects/004_Propofol/Modelling/neuron_simulations/code')
import niktools as nt
import numpy as np
from matplotlib import pyplot as plt
import LFPy
import fig1
from neuron import h
import mat4py as m4p


propofol = int(sys.argv[1])
neuronID = sys.argv[2]

################################################
################## Morphology ##################
################################################

filepath = 'E:/Research_Projects/004_Propofol/Modelling/neuron_simulations/data/dendrite_morphologies/prinicple_cells/'
file = filepath + neuronID

pts3d,connections,segs,morphData = fig1.morph2Segs(file)
dendIdcs = np.concatenate(segs)
p0 = np.mean(morphData[dendIdcs-1,2:5],axis=0)

l = list()
for i in range(len(segs)):
    s = np.array(morphData[segs[i]-1,3:6])
    l.append(np.sum(np.linalg.norm(np.diff(s,axis=0),2,1)))

L = np.sum(l)

Epas = -80 # a la Ahlfors and Wreh
Gpas = 1/30000
Ra = 150 

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


###################################################
################ Set up cell model ################
###################################################
tmax = 2000
tV = np.linspace(0,tmax,10000)

iSyn = {'tau1':1,
        'tau2':10*(1+propofol), # 10 ms
        'weight':0.1,#014*(1+propofol), # uS (1.4 nS) [Zsiros and Hestrin. J. Neurophysiol. 2005] (10-30 pS)
        'erev':-80,
        'lambda': 0.3*2*(1-0.3*propofol)}
eSyn = {'tau1':0.5,
        'tau2':3,
        'weight':0.1/100,#007, # uS (0.7 ns) [Zsiros and Hestrin. J. Neurophysiol. 2005]
        'erev':0,
        'lambda':0.3*0.5*(1-0.3*propofol)}


synType = 'Exp2Syn'
synParamsI = {'idx': 0,'e': iSyn['erev'],'syntype': synType,'tau1': iSyn['tau1'],'tau2': iSyn['tau2'],'weight': iSyn['weight'],'record_current': False}
synParamsE = {'idx': 0,'e': eSyn['erev'],'syntype': synType,'tau1': 1,'tau2': 3,'weight': eSyn['weight'],'record_current': False}

# synType = 'ExpSyn'
# synParamsI = {'idx': 0,'e': iSyn['erev'],'syntype': synType,'tau': iSyn['tau2'],'weight': iSyn['weight'],'record_current': False}
# synParamsE = {'idx': 0,'e': eSyn['erev'],'syntype': synType,'tau': eSyn['tau2'],'weight': eSyn['weight'],'record_current': False}


cell = LFPy.Cell(sectionList,v_init=Epas,celsius=37)
cell.tstop = tmax
syn = list()
synIds = list()
ei = list()
spike_times = list()

############### Excitatory synapses ###############
lambdaE = eSyn['lambda']  # rate per synapse
fE = lambda t,L: 1-0*.9*(np.sin(2*np.pi*t*5/1000)>0)
yE = np.cumsum(fE(tV,lambdaE))/np.sum(fE(tV,lambdaE))
nE = 100
tE = (500*np.ones(nE)).tolist()

synE = cell.get_rand_idx_area_norm(section='allsec', nidx=nE, z_min=575, z_max=1000000.0)

tsE = list()
for i in range(nE):
    tsE.append([])

for i in tE:
    j = np.random.randint(nE)
    tsE[j].append(i)

for i in range(nE):
    synParams = synParamsE.copy()
    synParams['idx'] = synE[i].tolist()
    syn.append(LFPy.Synapse(cell, **synParams))
    syn[-1].set_spike_times(np.array(tsE[i]))
    synIds.append(synE[i].tolist())
    ei.append('e')
    x.append(np.array([syn[-1].x,syn[-1].y,syn[-1].z]).tolist() )
    spike_times.append(tsE[i])


############### Inhibitory synapses ###############
lambdaI = iSyn['lambda']
fI = lambda t,L: 1-0*.9*(np.sin(2*np.pi*t*5/1000)>0)
yI = np.cumsum(fI(tV,lambdaI))/np.sum(fI(tV,lambdaI))
nI = 0
tI = (500*np.ones(nI)).tolist()

synI = cell.get_rand_idx_area_norm(section='allsec', nidx=nI, z_min=-400, z_max=0)


tsI = list()
for i in range(nI):
    tsI.append([])

for i in tI:
    j = np.random.randint(nI)
    tsI[j].append(i)

for i in range(nI):
    synParams = synParamsI.copy()
    synParams['idx'] = synI[i].tolist()
    syn.append(LFPy.Synapse(cell, **synParams))
    syn[-1].set_spike_times(np.array(tsI[i]))
    synIds.append(synI[i].tolist())
    ei.append('i')
    x.append(np.array([syn[-1].x,syn[-1].y,syn[-1].z]).tolist() )
    spike_times.append(tsI[i])

############### Simulate model ###############
cdm = LFPy.CurrentDipoleMoment(cell=cell)
cell.simulate(probes=[cdm],rec_vmem= True)

################ Save output ################
V = cell.vmem[0,:]
Q = cdm.data
Qsave = [Q[0,:].T.tolist(),Q[2,:].T.tolist(),(-Q[1,:].T).tolist()]


synData = {'e_properties': eSyn,
            'i_properties': iSyn,
            'ids':synIds,
            'ei':ei,
            'locations':x,
            'times':spike_times}

morphology = {'segments':[x.tolist() for x in segs],
            'coordinates':morphData.tolist(),
            'Vx':cell.x.tolist(),
            'Vy':cell.y.tolist(),
            'Vz':cell.z.tolist(),
            'fileName':file}

data = {'Epas':Epas,
    'Ra':Ra,
    't':cell.tvec.tolist(),
    'synapses':synData,
    'morphology':morphology,
    'V':V.tolist(),
    'Q':Qsave}


savePath = 'E:/Research_Projects/004_Propofol/Modelling/neuron_simulations/data/simulations/temp/'
# savePath = 'E:/Research_Projects/004_Propofol/Modelling/neuron_simulations/data/simulations/before_after_propofol/2022_04_27 parameters/'
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

# nt.bell()

print(savedName)