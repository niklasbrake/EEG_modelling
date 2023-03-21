## ------------ UPDATE 2022-05-19 ------------
## 1.   Turned firing rates of presyanptic neurons to be log normal with median
##      equal to prevous rates and sigma = 0.4
##
## 2.   Removed excitatory syanpses on the soma
##
## 3.   Added correlation in inputs between each N neurons. Synapse placements now 
##      the SAME on each neuron repetition

# savePath = 'E:/Research_Projects/004_Propofol/Modelling/neuron_simulations/data/simulations/examples/'
savePath = 'E:/Research_Projects/004_Propofol/Modelling/neuron_simulations/data/simulations/cortical_column/baseline_correlation/';

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

load_mechanisms('E:/Research_Projects/004_Propofol/Modelling/neuron_simulations/code/_mod')

propofol = int(sys.argv[1])
neuronID = sys.argv[2]
nReps = int(sys.argv[3])
corr = float(sys.argv[4])

print(corr)

# gamShape = 512 # Poisson
gamShape = 1 # Poisson

################################################
################## Morphology ##################
################################################

# filepath = 'E:/Research_Projects/004_Propofol/Modelling/neuron_simulations/data/dendrite_morphologies/prinicple_cells/'
filepath = 'E:/Research_Projects/004_Propofol/Modelling/neuron_simulations/data/cortical_column_Hagen/swc/'
file = filepath + neuronID

pts3d,connections,segs,morphData = getMorphoSegments.morph2Segs(file)
dendIdcs = np.concatenate(segs)
p0 = np.mean(morphData[dendIdcs-1,2:5],axis=0)

l = list()
for i in range(len(segs)):
    s = np.array(morphData[segs[i]-1,3:6])
    l.append(np.sum(np.linalg.norm(np.diff(s,axis=0),2,1)))

L = np.sum(l)

Ra = 100 
eLeak = -63 # Masoli et al
gLeak = 0.0011

# Add soma
soma = h.Section(name='soma')

soma.insert('Leak')
soma.e_Leak = eLeak
soma.gmax_Leak = gLeak

# soma.insert('inaT')
# soma.shift_inaT = -10
# soma.vtraub_inaT = -70
# soma.gnabar_inaT = 0.12*5
# soma.ena = 90

# soma.insert('ikdT')
# soma.gkbar_ikdT = 0.1*5
# soma.ek = -80

# soma.insert('imZ')
# soma.gkbar_imZ = 5*1e-5

soma.push()
h.pt3dclear()
i = 0
while morphData[i,1]==1:
    h.pt3dadd(morphData[i,2], morphData[i,3], morphData[i,4], morphData[i,5])
    i+=1

sectionList = h.SectionList()
sectionList.append()
h.pop_section()

# soma.nseg = 5

# Add 3D points for dendritic tree
dend = list()
x = list()
for i in range(len(segs)):
    dend.append(h.Section(name='dend'+str(i)))
    dend[-1].push()
    h.pt3dclear()
    for j in segs[i]:
        idx = np.argwhere(morphData[:,0]==j)[0][0]
        h.pt3dadd(morphData[idx,2],morphData[idx,3],morphData[idx,4],morphData[idx,5])
    # h.pt3dadd(pts3d[i,0], pts3d[i,1], pts3d[i,2], pts3d[i,3])
    # h.pt3dadd(pts3d[i,4], pts3d[i,5], pts3d[i,6], pts3d[i,7])
    dend[-1].insert('Leak')
    dend[-1].e_Leak = eLeak
    dend[-1].gmax_Leak = gLeak
    # dend[-1].nseg = 3
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


####################################################
################ Syanpse parameters ################
####################################################
eSyn = {'tau1':0.3,
        'tau2':1.8,
        'weight':0.0007, # uS (0.7 ns) [Gidon et al., 2020]
        'erev':0,
        'lambda': 0.3*1*(1-0.5*propofol)}
iSyn = {'tau1':1,
        'tau2':23*(1+propofol), # 10 ms
        'tau1_proximal':0.5, # 10 ms
        'tau2_proximal':5*(1+propofol), # 10 ms
        'weight':0.0014*(1+propofol), # uS (1.4 nS) [Gidon et al., 2020] (10-30 pS)
        'erev':-75,
        'lambda': 0.3*5*(1-0*0.3*propofol)}

synType = 'Exp2Syn'
synParamsE = {'idx': 0,
                'e': eSyn['erev'],
                'syntype': synType,
                'tau1': eSyn['tau1'],
                'tau2': eSyn['tau2'],
                'weight': eSyn['weight'],
                'record_current': False}

synParamsI = {'idx': 0,
                'e': iSyn['erev'],
                'syntype': synType,
                'tau1': iSyn['tau1'],
                'tau2': iSyn['tau2'],
                'weight': iSyn['weight'],
                'record_current': False}

tmax = 2000



############### Excitatory synapses ##############
# Get location of each excitatory syanpse
T = eSyn['lambda']*tmax/1000 # total number events
maxN = int(poisson.ppf(0.9999,T))
nE = np.random.poisson(L) # total number of syanpses (1/um)

cell = LFPy.Cell(sectionList,v_init=eLeak,celsius=37)
synE = cell.get_rand_idx_area_norm(section='allsec', nidx=nE, z_min=- 1000000.0, z_max=1000000.0)

def getSpikeTimes(gamShape,lam,maxN):
    t0 = 1000*np.cumsum(np.random.gamma(gamShape,1/gamShape/lam,maxN))
    offset = np.random.rand(len(t0))*1000/lam
    ts = [t-b for t,b in zip(t0,offset) if t-b<=tmax and t-b >= 0] 
    return ts

# Remove somatic Ex connections
idcs = np.argwhere(cell.get_idx_name(synE)[:,1]=='soma')
while len(idcs)>0:
    synE[idcs] = cell.get_rand_idx_area_norm(section='allsec', nidx=len(idcs), z_min=- 1000000.0, z_max=1000000.0)[:,None]
    idcs = np.argwhere(cell.get_idx_name(synE)[:,1]=='soma')



tsE_common = getSpikeTimes(gamShape,corr*eSyn['lambda'],maxN)

# lognormal distribution of presynaptic firing rates
lambdaSynE = list()
for i in range(nE*nReps):
    lambdaSynE.append(np.random.lognormal(np.log(eSyn['lambda']),0.4,1)[0])
    tsE = getSpikeTimes(gamShape,(1-corr)*lambdaSynE[-1],maxN)
    popDriveT = [x for x in tsE_common]
    tsE.extend(popDriveT)

# Generate excitatory spike times for each presynaptic neuron
tsE = getSpikeTimes(gamShape,(1-corr)*lambdaSynE[i],maxN)
preESynTimes = list()
for i in range(nE):
    temp = list()
    tsE_common = getSpikeTimes(gamShape,corr*lambdaSynE[i],maxN)
    for rep in range(nReps):
        tsE.extend(tsE_common)
        temp.append(tsE)
    preESynTimes.append(temp)
ESynTimes = list(map(list,zip(*preESynTimes))) # transpose list
###################################################

############### Inhibitory synapses ###############
# Get location of each inhibitory syanpse
nI = np.random.poisson(0.15*L)  # total number of syanpses (0.15/um)
synI = cell.get_rand_idx_area_norm(section='allsec', nidx=nI, z_min=-1000, z_max=10000)

T = iSyn['lambda']*tmax/1000 # total number events
maxN = int(poisson.ppf(0.9999,T))

# lognormal distribution of presynaptic firing rates
lambdaSynI = list()
for i in range(nE):
    lambdaSynI.append(np.random.lognormal(np.log(iSyn['lambda']),0.4,1)[0])

# Generate inhibitory spike times for each presynaptic neuron
preISynTimes = list()
tsI_common = getSpikeTimes(gamShape,corr*lambdaSynI[i],maxN)
for i in range(nI):
    temp = list()
    for rep in range(nReps):
        tsI = getSpikeTimes(gamShape,(1-corr)*lambdaSynI[i],maxN)
        tsI.extend(tsI_common)
        temp.append(tsI)
    preISynTimes.append(temp)
ISynTimes = list(map(list,zip(*preISynTimes))) # transpose list
###################################################

# Loop through each repeition, with different spike times
for rep in range(nReps):
    print(rep)
    # np.random.seed()
    cell = LFPy.Cell(sectionList,v_init=eLeak,celsius=37)
    cell.tstop = tmax
    syn = list()
    synIds = list()
    ei = list()
    spike_times = list()
    x = list()

    # set spike times for every synapse
    for i in range(nE):
        synParams = synParamsE.copy()
        synParams['idx'] = synE[i].tolist()
        syn.append(LFPy.Synapse(cell, **synParams))
        synIds.append(synE[i].tolist())
        ei.append('e')
        x.append(np.array([syn[-1].x,syn[-1].y,syn[-1].z]).tolist() )

        syn[-1].set_spike_times(np.array(ESynTimes[rep][i]))
        spike_times.append(ESynTimes[rep][i])

    for i in range(nI):
        synParams = synParamsI.copy()
        synParams['idx'] = synI[i].tolist()
        syn.append(LFPy.Synapse(cell, **synParams))
        synIds.append(synI[i].tolist())
        ei.append('i')
        x.append(np.array([syn[-1].x,syn[-1].y,syn[-1].z]).tolist())

        syn[-1].set_spike_times(np.array(ISynTimes[rep][i]))
        spike_times.append(ISynTimes[rep][i])

    ############### Simulate model ###############
    cdm = LFPy.CurrentDipoleMoment(cell=cell)
    cell.simulate(probes=[cdm],rec_somav=True,rec_vmem=True)

    ################ Save output ################
    V = cell.vmem
    V_soma = cell.somav
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

    data = {'eLeak':eLeak,
        'Ra':Ra,
        't':cell.tvec.tolist(),
        'synapses':synData,
        'morphology':morphology,
        'V':V.tolist(),
        'V_soma':V_soma.tolist(),
        'Q':Qsave}

    
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


    print(savedName)

nt.bell()