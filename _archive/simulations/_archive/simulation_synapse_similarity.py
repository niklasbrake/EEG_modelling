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
    preEI = list()
    preSpikes = list()
    with open(spikingFile, newline='') as csvfile:
         spikeTimes = csv.reader(csvfile, delimiter=' ', quotechar='|')
         for row in spikeTimes:
            temp = row[0].split(',')
            if(temp[0][0]=='e'):
                preEI.append(0)
                if(temp[1:]==['']):
                    preSpikes.append(list())
                else:
                    preSpikes.append([float(x) for x in temp[1:]])
            elif(temp[0][0]=='i'):
                preEI.append(1)
                if(temp[1:]==['']):
                    preSpikes.append(list())
                else:
                    preSpikes.append([float(x) for x in temp[1:]])
    return preEI,preSpikes

def load_syanpse_locations(synapseFile):
    """ Load synapse locations from the file specified
    by the variable synapseFile """
    synSeg = list()
    synPre = list()
    with open(synapseFile, newline='') as csvfile:
         syns = csv.reader(csvfile, delimiter=' ', quotechar='|')
         for row in syns:
            temp = row[0].split(',')
            synSeg.append(int(temp[0])-1)
            synPre.append(int(temp[1])-1)
    return synSeg,synPre

def save_simulation(data,sim_path,nrn_name,propofol):
    """ Save simulation results as a matlab file """
    if(propofol):
        condition = '_propofol_simulation'
    else:
        condition = '_simulation'
    number = 1
    savedName = sim_path + '/LFPy/'  + nrn_name + condition + f"{number:02d}" + '.mat'
    while(os.path.exists(savedName)):
        number += 1
        savedName = sim_path + '/LFPy/'  + nrn_name + condition + f"{number:02d}" + '.mat'
    m4p.savemat(savedName, data)
    print(savedName)

def get_area_distribution(neuronID):
    AXIAL_RESISTANCE        = 100
    EREV_LEAK               = -63       # Masoli et al
    G_LEAK                  = 0.0011

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
                    'tau2': 23,
                    'weight': 0.0014, # uS
                    'record_current': False}

    filepath = ('E:/Research_Projects/004_Propofol/Modelling/'
        'neuron_simulations/data/cortical_column_Hagen/swc/')
    file = filepath + neuronID + '.swc'
    pts3d,connections,segs,morphData = getMorphoSegments.morph2Segs(file)

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

    cell = LFPy.Cell(sectionList,v_init=EREV_LEAK,celsius=37)
    savedName = 'E:/Research_Projects/004_Propofol/Modelling/neuron_simulations/data/cortical_column_Hagen/segment_area_distributions/' + neuronID + '_area.mat'
    data = {'area':cell.area.tolist(),
            'x':cell.x.tolist(),
            'y':cell.y.tolist(),
            'z':cell.z.tolist()}
    m4p.savemat(savedName, data)

def addsyns(synapseFile,preEI,preSpikes,propofol):
    SYANPSE_FAILURE_RATE    = 0 # 0.7;      # Synapse failure rate
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

    synSeg,synPre = load_syanpse_locations(synapseFile)
    synTimes = list()
    synParams = list()
    # Get excitatory syanpse locations
    for i in range(len(synSeg)):
        ts = preSpikes[synPre[i]]
        spikeSelection = np.random.binomial(len(ts),1-SYANPSE_FAILURE_RATE)
        ts = np.random.choice(ts,spikeSelection,replace=False).tolist()
        if(preEI[synPre[i]]):
            params = IN_PARAMS.copy()
        else:
            params = EX_PARAMS.copy()
        params['idx'] = synSeg[i]
        synParams.append(params)
        synTimes.append(ts)

    return synParams,synTimes

def main(propofol,neuronID,nReps,sim_path,activeSoma=False,T_MAX=2000):
    AXIAL_RESISTANCE        = 100
    EREV_LEAK               = -63       # Masoli et al
    G_LEAK                  = 0.0011

    filepath = ('E:/Research_Projects/004_Propofol/Modelling/'
        'neuron_simulations/data/cortical_column_Hagen/swc/')
    file = neuronID + '.swc'
    pts3d,connections,segs,morphData = getMorphoSegments.morph2Segs(filepath+file)
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

    spikingFile = sim_path + '/network_spikes.csv'
    preEI,preSpikes = load_spike_times(spikingFile)

    # Simulate nReps iterations of the model
    for rep in range(nReps):
        cell = LFPy.Cell(sectionList,v_init=-65,celsius=37)
        cell.tstop = T_MAX

        synapseFile = sim_path + '/network/' + neuronID + '_N' + str(rep+1) + '.csv'
        synParams,synTimes = addsyns(synapseFile,preEI,preSpikes,propofol)
        syn = list()
        for i in range(len(synParams)):
            syn.append(LFPy.Synapse(cell, **synParams[i]))
            syn[-1].set_spike_times(np.array(synTimes[i]))

        cdm = LFPy.CurrentDipoleMoment(cell=cell)
        cell.simulate(probes=[cdm],rec_somav=True,rec_vmem=False)

        Qx = cdm.data[0,:].T.tolist()
        Qy = cdm.data[1,:].T.tolist()
        Qz = cdm.data[2,:].T.tolist()

        data = {'t':cell.tvec.tolist(),
            'V_soma':cell.somav.tolist(),
            'Q':[Qx,Qy,Qz],
            'synapseFile':synapseFile}

        save_simulation(data,sim_path,neuronID,propofol)

    nt.bell()
    return data

if __name__ == "__main__":
    n = len(sys.argv)
    propofol = int(sys.argv[1])
    neuronID = sys.argv[2]
    if(n==3):
        get_area_distribution(neuronID)
    else:
        nReps = int(sys.argv[3])
        sim_path = sys.argv[4]
        if(n==6):
            T_MAX = int(sys.argv[5])
        else:
            T_MAX = 2000
        main(propofol,neuronID,nReps,sim_path,False,T_MAX)
