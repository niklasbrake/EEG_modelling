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
sys.path.append('C:/Users/brake/Documents/GitHub/aperiodic_EEG_modelling/simulations/functions')
import niktools as nt
import numpy as np
import LFPy
import getMorphoSegments
from cortical_neuron import cortical_neuron as init_neuron
import csv
import mat4py as m4p
from neuron import h, load_mechanisms
from matplotlib import pyplot as plt
from scipy.stats import poisson
from math import ceil

def import_postsyanptic_network(networkPath):
    file = networkPath + '/mTypes.txt'
    mFiles = list()
    with open(file, newline='') as csvfile:
         syns = csv.reader(csvfile, delimiter='\n', quotechar='|')
         for row in syns:
            mFiles.append(row[0])
    synDir = networkPath + '/connections'
    synapseFiles = [synDir + '/' + x for x in os.listdir(synDir)]
    return mFiles,synapseFiles

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
    # SYANPSE_FAILURE_RATE    = 0.7;      # Synapse failure rate

    # tau = int(propofol/100)
    # remainder = propofol-tau*100
    # tau_change = int(remainder/10)
    # weight_change = remainder-tau_change*10
    EX_PARAMS = {'idx': 0,
                    'e': 0,
                    'syntype': 'Exp2Syn',
                    'tau1': 0.3,
                    'tau2': 1.8,
                    'weight': 0.0007, # uS
                    'record_current': False}
    IN_PARAMS = {'idx': 0,
                    # 'e': -75,
                    'e':-80,
                    'syntype': 'Exp2Syn',
                    'tau1': 2,
                    'tau2': 23*(1+propofol),
                    # 'tau2': tau*tau_change,
                    # 'tau2': 16*(1+propofol),
                    'weight': 0.0014,#*weight_change, # uS
                    'record_current': False}
    synSeg,synPre = load_syanpse_locations(synapseFile)
    synTimes = list()
    synParams = list()
    # Get excitatory syanpse locations
    for i in range(len(synSeg)):
        ts = preSpikes[synPre[i]]
        # ts = [x+1000 for x in ts]
        # ts = [x for x in ts if x>1e3]
        # spikeSelection = np.random.binomial(len(ts),1-SYANPSE_FAILURE_RATE)
        # ts = np.random.choice(ts,spikeSelection,replace=False).tolist()
        if(preEI[synPre[i]]):
            params = IN_PARAMS.copy()
        else:
            params = EX_PARAMS.copy()
        params['idx'] = synSeg[i]
        synParams.append(params)
        synTimes.append(ts)

    return synParams,synTimes

def main(mFile,synapseFiles,preEI,preSpikes,savePath,T_MAX=100,activeSoma=False,propofol=0):

    # Initialize neuron morphology
    nrnM = init_neuron(mFile,activeSoma)

    # Simulate neurons with given morphology
    for file in synapseFiles:
        cell = LFPy.Cell(nrnM.sectionList,v_init=-65,celsius=37)
        cell.tstop = T_MAX

        # Add synapses to neuron model
        synParams,synTimes = addsyns(file,preEI,preSpikes,propofol)
        syn = list()
        for i in range(len(synParams)):
            syn.append(LFPy.Synapse(cell, **synParams[i]))
            syn[-1].set_spike_times(np.array(synTimes[i]))

        # Simulate
        cdm = LFPy.CurrentDipoleMoment(cell=cell)
        cell.simulate(probes=[cdm],rec_somav=True,rec_vmem=False)

        #Save results
        nrnID = str.split(file,'/')[-1][:-4]
        saveFile = savePath+'/LFPy/'+nrnID +'.csv'
        t = cell.tvec.reshape([-1,1])
        v = cell.somav.reshape([-1,1])
        data = np.concatenate((t,cdm.data.T,v),axis=1)
        np.savetxt(saveFile, data, delimiter=",", header="time,Qx,Qy,Qz,V_soma", fmt="%f,%f,%f,%f,%f")
        print(saveFile)
        cell.strip_hoc_objects()


if __name__ == "__main__":
    pars = sys.argv
    n = len(pars)
    networkPath = pars[1]
    spikingFile = pars[2]
    savePath = pars[3]
    T_MAX = int(pars[4])
    activeSoma = (pars[5]=="true")
    propofol = float(pars[6])

    preEI,preSpikes = load_spike_times(spikingFile)
    mFiles,synapseFiles = import_postsyanptic_network(networkPath)

    mTypes = dict()
    for i,m in enumerate(mFiles):
        if m not in mTypes:
            mTypes[m] = [synapseFiles[i]]
        else:
            mTypes[m].append(synapseFiles[i])

    for m in mTypes.keys():
        main(m,mTypes[m],preEI,preSpikes,savePath,T_MAX,activeSoma,propofol)