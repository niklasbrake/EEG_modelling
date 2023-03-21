"""
------------ UPDATE 2022-05-19 ------------
Removed excitatory syanpses on the soma

------------ UPDATE 2022-05-30 ------------
Presynaptic spike times are now an input, and synapse activaiton is
drawn from these spike times. 70% of spikes fail to activate synapse

------------ UPDATE 2022-06-01 ------------
Three seperate populations. One for apical dendrites, another for basal,
and a third for soma. Neurons that doesn't have apical dendrites will get
input from both apical and basal presynaptic populations.

------------ UPDATE 2022-06-17 ------------
Different time scales for inhibitory syanpses onto soma and basal
dendrites versus apical dendrites, a la Gidon et al., Science 2020.
Also decreased the untiary conductance from 1.4 -> 0.5 following
Gidon et al.
"""

import sys
import os
sys.path.append('C:/Users/brake/Documents/pyNB')
sys.path.append('E:/Research_Projects/004_Propofol/Modelling/neuron_simulations')
sys.path.append('E:/Research_Projects/004_Propofol/Modelling/neuron_simulations/code')
import niktools as nt
import numpy as np
import LFPy
import getMorphoSegments
import csv
import mat4py as m4p
from neuron import h, load_mechanisms
from matplotlib import pyplot as plt
from scipy.stats import poisson
from math import ceil

def load_spike_times(spikingFile):
    """ Load spike times from the file specified
    by the variable spikingFile """
    exApical = list()
    exBasal = list()
    inSoma = list()
    inApical = list()
    inBasal = list()
    with open(spikingFile, newline='') as csvfile:
         spikeTimes = csv.reader(csvfile, delimiter=' ', quotechar='|')
         for row in spikeTimes:
            temp = row[0].split(',')
            if(temp[0]=='eA'):
                if(temp[1:]==['']):
                    exApical.append(list())
                else:
                    exApical.append([float(x) for x in temp[1:]])
            elif(temp[0]=='eB'):
                if(temp[1:]==['']):
                    exBasal.append(list())
                else:
                    exBasal.append([float(x) for x in temp[1:]])
            elif(temp[0]=='iA'):
                if(temp[1:]==['']):
                    inApical.append(list())
                else:
                    inApical.append([float(x) for x in temp[1:]])
            elif(temp[0]=='iB'):
                if(temp[1:]==['']):
                    inBasal.append(list())
                else:
                    inBasal.append([float(x) for x in temp[1:]])
            elif(temp[0]=='iS'):
                if(temp[1:]==['']):
                    inSoma.append(list())
                else:
                    inSoma.append([float(x) for x in temp[1:]])

    return exApical,exBasal,inSoma,inApical,inBasal

def save_simulation(data,file,propofol):
    """ Save simulation results as a matlab file """
    nrn_name = os.path.split(file)[-1][:-4]
    if(propofol):
        condition = '_propofol_simulation'
    else:
        condition = '_simulation'
    number = 1
    savedName = savePath  + '/' + nrn_name + condition + f"{number:02d}" + '.mat'
    while(os.path.exists(savedName)):
        number += 1
        savedName = savePath  + '/' + nrn_name + condition + f"{number:02d}" + '.mat'
    m4p.savemat(savedName, data)
    print(savedName)

def main(propofol,neuronID,nReps,spikingFile,savePath,activeSoma=False,T_MAX=2000):
    """
    propofol = 0
    neuronID = 'L23E_oi24rpy1'
    savePath = 'E:/Research_Projects/004_Propofol/Modelling/neuron_simulations/data/simulations/cortical_column/oscillations/20220615/dipoles'
    spikingFile = 'E:/Research_Projects/004_Propofol/Modelling/neuron_simulations/data/simulations/cortical_column/oscillations/20220615/spikeTimes.csv'
    nReps = 1
    activeSoma = True
    """
    AXIAL_RESISTANCE        = 100
    EREV_LEAK               = -63       # Masoli et al
    G_LEAK                  = 0.0011
    SYANPSE_FAILURE_RATE    = 0 # 0.7;      # Synapse failure rate
    SYNS_IN_CON             = 1         # or 13.9 from Markram et al.
    SYNS_EX_CON             = 1         # or 3.6 Markram et al.

    EX_PARAMS = {'idx': 0,
                    'e': 0,
                    'syntype': 'Exp2Syn',
                    'tau1': 0.3,
                    'tau2': 1.8,
                    'weight': 0.0007, # uS
                    'record_current': False}

    IN_PARAMS = {'idx': 0,
                    'e': -75,
                    'syntype': 'Exp2Syn',
                    'tau1': 2,
                    'tau2': 23*(1+propofol),
                    'weight': 0.0014*(1+propofol), # uS
                    'record_current': False}

    filepath = ('E:/Research_Projects/004_Propofol/Modelling/'
        'neuron_simulations/data/cortical_column_Hagen/swc/')
    file = filepath + neuronID + '.swc'
    pts3d,connections,segs,morphData = getMorphoSegments.morph2Segs(file)

    # get total length of dendrites
    l = list()
    for i in range(len(segs)):
        s = np.array(morphData[segs[i]-1,3:6])
        l.append(np.sum(np.linalg.norm(np.diff(s,axis=0),2,1)))

    L_arbour = np.sum(l)

    load_mechanisms('E:/Research_Projects/004_Propofol/Modelling/neuron_simulations/code/_mod')

    # List to record the type of neuron segment:
    # soma = 0, basal dendrite = 3, apical dendrite = 4
    secType = dict()

    # Generate NEURON list of sections
    sectionList = h.SectionList()

    # Add soma to NEURON model
    soma = h.Section(name='soma')
    secType['soma'] = 1
    soma.insert('Leak')
    soma.e_Leak = EREV_LEAK
    soma.gmax_Leak = G_LEAK
    if activeSoma:
        soma.insert('inaT')
        soma.shift_inaT = 1
        soma.vtraub_inaT = -74
        soma.gnabar_inaT = 0.12*5
        soma.ena = 60
        soma.insert('ikdT')
        soma.gkbar_ikdT = 0.1*5
        soma.ek = -80
        # soma.insert('imZ')
        # soma.gkbar_imZ = 5*1e-5

    soma.push()
    h.pt3dclear()
    i = 0
    while morphData[i, 1] == 1:
        h.pt3dadd(morphData[i,2], morphData[i,3], morphData[i,4], morphData[i,5])
        i += 1

    sectionList.append()
    h.pop_section()

    # Add dendrite segments to NEURON model
    dend = list()
    for i in range(len(segs)):
        secName = 'dend'+str(i)
        dend.append(h.Section(name=secName))
        dend[-1].push()
        h.pt3dclear()
        dendType = list()
        for j in segs[i]:
            idx = np.argwhere(morphData[:,0]==j)[0][0]
            h.pt3dadd(morphData[idx,2],morphData[idx,3],morphData[idx,4],morphData[idx,5])
            dendType.append(morphData[idx,1])
        if(len(dendType)==2):
            secType[secName] = dendType[-1]
        else:
            secType[secName] = max(set(dendType), key=dendType.count)
        dend[-1].insert('Leak')
        dend[-1].e_Leak = EREV_LEAK
        dend[-1].gmax_Leak = G_LEAK
        dend[-1].Ra = AXIAL_RESISTANCE
        sectionList.append()
        h.pop_section()

    # Connect all compartments together
    for i in range(len(connections)):
        i1 = connections[i][0]
        i0 = connections[i][1]
        if(i1==-1):
            dend[i0].connect(soma,1,0)
        else:
            dend[i0].connect(dend[i1](1))

    exApical,exBasal,inSoma,inApical,inBasal = load_spike_times(spikingFile)

    # Simulate nReps iterations of the model
    for rep in range(nReps):
        cell = LFPy.Cell(sectionList,v_init=EREV_LEAK,celsius=37)
        cell.tstop = T_MAX
        syn = list()
        synID = list()
        synType = list()
        synTimes = list()
        synXYZ = list()

        # Get excitatory syanpse locations
        nE = np.random.poisson(L_arbour)
        synE = cell.get_rand_idx_area_norm(section='allsec', nidx=nE, z_min=- 1000000.0, z_max=1000000.0)

        # Remove somatic Ex connections
        idcs = np.argwhere(cell.get_idx_name(synE)[:,1]=='soma')
        while len(idcs)>0:
            synE[idcs] = cell.get_rand_idx_area_norm(section='allsec', nidx=len(idcs), z_min=- 1000000.0, z_max=1000000.0)[:,None]
            idcs = np.argwhere(cell.get_idx_name(synE)[:,1]=='soma')

        eTargets = [secType[cell.get_idx_name(x)[1]] for x in synE]
        nApical = sum([x==4 for x in eTargets])
        nBasal = sum([x==3 for x in eTargets])
        nSomatic = sum([x==1 for x in eTargets])
        if(nSomatic>0):
            print('Error: excitatory somatic input');
        # Randomly draw presynaptic neurons
        eApicalSpikeTrains = np.random.choice(exApical,int(nApical/SYNS_EX_CON),replace=False)
        if(nApical==0): # if no apical dendrites
            eBasalSpikeTrains = np.random.choice(exBasal+exApical,int(nBasal/SYNS_EX_CON),replace=False)
        else:
            eBasalSpikeTrains = np.random.choice(exBasal,int(nBasal/SYNS_EX_CON),replace=False)

        for i in range(nE):
            # Add synapse to NEURON model and save location
            synParams = EX_PARAMS.copy()
            synParams['idx'] = synE[i].tolist()
            syn.append(LFPy.Synapse(cell, **synParams))
            synID.append(synE[i].tolist())
            synType.append('e')
            synXYZ.append(np.array([syn[-1].x,syn[-1].y,syn[-1].z]).tolist())

            # Get spike train from random presynaptic neurons
            dendType = secType[cell.get_idx_name(synE[i])[1]] # Type of dendrite
            if(dendType==0): # Soma
                ts = np.random.choice(eBasalSpikeTrains) # All somatic esyns should have been removed
            elif(dendType==3): # Basal dendrites
                ts = np.random.choice(eBasalSpikeTrains)
            elif(dendType==4): # Apical dendrites
                ts = np.random.choice(eApicalSpikeTrains)
            synapseTimes = np.random.binomial(len(ts),1-SYANPSE_FAILURE_RATE)
            ts = np.random.choice(ts,synapseTimes,replace=False).tolist()
            syn[-1].set_spike_times(np.array(ts))
            synTimes.append(ts)

        # Get inhibitory syanpse locations
        nI = np.random.poisson(0.15*L_arbour)
        synI = cell.get_rand_idx_area_norm(section='allsec', nidx=nI, z_min=-1000, z_max=10000)

        # Randomly draw presynaptic neurons
        iTargets = [secType[cell.get_idx_name(x)[1]] for x in synI]
        nApical = sum([x==4 for x in iTargets])
        nBasal = sum([x==3 for x in iTargets])
        nSomatic = sum([x==1 for x in iTargets])
        iSomaticSpikeTrains = np.random.choice(inSoma,int(nSomatic/SYNS_IN_CON),replace=False)
        iApicalSpikeTrains = np.random.choice(inApical,int(nApical/SYNS_IN_CON),replace=False)
        if(nApical==0): # if no apical dendrites
            iBasalSpikeTrains = np.random.choice(inBasal+inApical,int(nBasal/SYNS_IN_CON),replace=False)
        else:
            iBasalSpikeTrains = np.random.choice(inBasal,int(nBasal/SYNS_IN_CON),replace=False)

        for i in range(nI):
            # Get spike train from random presynaptic neurons
            dendType = secType[cell.get_idx_name(synE[i])[1]] # Type of dendrite
            if(dendType==0): # Soma
                synParams = IN_PARAMS.copy()
                ts = np.random.choice(iSomaticSpikeTrains)
            elif(dendType==3): # Basal dendrites
                synParams = IN_PARAMS.copy()
                ts = np.random.choice(iBasalSpikeTrains)
            elif(dendType==4): # Apical dendrites
                synParams = IN_PARAMS.copy()
                ts = np.random.choice(iApicalSpikeTrains)
            synapseTimes = np.random.binomial(len(ts),1-SYANPSE_FAILURE_RATE)
            ts = np.random.choice(ts,synapseTimes,replace=False).tolist()
            # Add synapse to NEURON model and save location
            synParams['idx'] = synI[i].tolist()
            syn.append(LFPy.Synapse(cell, **synParams))
            syn[-1].set_spike_times(np.array(ts))
            # Save parameters for syanpse
            synID.append(synI[i].tolist())
            synType.append('i')
            synXYZ.append(np.array([syn[-1].x,syn[-1].y,syn[-1].z]).tolist())
            synTimes.append(ts)

        cdm = LFPy.CurrentDipoleMoment(cell=cell)
        cell.simulate(probes=[cdm],rec_somav=True,rec_vmem=False)

        Qx = cdm.data[0,:].T.tolist()
        Qy = cdm.data[1,:].T.tolist()
        Qz = cdm.data[2,:].T.tolist()

        morphology = {'segments':[x.tolist() for x in segs],
                    'coordinates':morphData.tolist(),
                    'Vx':cell.x.tolist(),
                    'Vy':cell.y.tolist(),
                    'Vz':cell.z.tolist(),
                    'fileName':file}

        synData = {'e_properties': EX_PARAMS,
            'i_properties': IN_PARAMS,
            'ids':synID,
            'synType':synType,
            'locations':synXYZ,
            'times':synTimes}

        data = {'EREV_LEAK':EREV_LEAK,
            'AXIAL_RESISTANCE':AXIAL_RESISTANCE,
            't':cell.tvec.tolist(),
            'synapses':synData,
            'morphology':morphology,
            'V':list(),
            'V_soma':cell.somav.tolist(),
            'Q':[Qx,Qy,Qz]}

        save_simulation(data,file,propofol)

    nt.bell()
    return data

if __name__ == "__main__":
    n = len(sys.argv)
    propofol = int(sys.argv[1])
    neuronID = sys.argv[2]
    nReps = int(sys.argv[3])
    spikingFile = sys.argv[4]
    savePath = sys.argv[5]
    if(n==7):
        T_MAX = int(sys.argv[6])
    else:
        T_MAX = 2000
    main(propofol,neuronID,nReps,spikingFile,savePath,False,T_MAX)
