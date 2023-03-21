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

def main(mFile,spikingFile,synapseFile,T_MAX=2000,activeSoma=False,propofol=0):
    AXIAL_RESISTANCE        = 100
    EREV_LEAK               = -63       # Masoli et al
    G_LEAK                  = 0.0011

    pts3d,connections,segs,morphData = getMorphoSegments.morph2Segs(mFile)
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

    preEI,preSpikes = load_spike_times(spikingFile)

    # Simulate nReps iterations of the model
    cell = LFPy.Cell(sectionList,v_init=-65,celsius=37)
    cell.tstop = T_MAX

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

    return data


def mFiles,syanpseFiles = import_postsyanptic_network(networkPath):
    file = networkPath + 'mTypes.txt'
    mFiles = list()
    with open(file, newline='') as csvfile:
         syns = csv.reader(csvfile, delimiter=' ', quotechar='|')
         for row in syns:
            mFiles.append(row)
    synDir = networkPath + '/connections'
    synapseFiles = [synDir + '/' + x for x in os.listdir(synDir)]

if __name__ == "__main__":
    pars = sys.argv
    n = len(pars)
    mFile = pars[1]
    spikingFile = pars[2]
    synapseFile = pars[3]
    savedName = pars[4]

    if n>4:
        data = main(mFile,spikingFile,synapseFile)
    if n>5:
        data = main(mFile,spikingFile,synapseFile,int(pars[5]))
    if n>6:
        data = main(mFile,spikingFile,synapseFile,int(pars[5]),pars[6])
    if n>7:
        data = main(mFile,spikingFile,synapseFile,int(pars[5]),pars[6],pars[7])
    m4p.savemat(savedName, data)
