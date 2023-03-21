import sys
import os
import numpy as np
import LFPy
import json
from cortical_neuron import cortical_neuron as init_neuron
import csv
from neuron import h, load_mechanisms
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

def addsyns(nrnID,meta_data,preEI,preSpikes,propofol):
    SYANPSE_FAILURE_RATE = 0.7
    if(propofol>100):
        tau = int(propofol/100)
        remainder = propofol-tau*100
        tau_change = int(remainder/10)
        weight_change = remainder-tau_change*10
    else:
        tau = 23
        tau_change = 1+2*propofol
        weight_change = 1+0*propofol

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
                    'tau1': 1,
                    # 'tau2': 13*(1+propofol*2),
                    'tau2': tau*tau_change,
                    # 'weight': 0.0014*(1+propofol*0), # uS
                    'weight': 0.0007*weight_change, # uS
                    'record_current': False}

    synSeg = meta_data['synSeg'][nrnID]
    synPre = meta_data['synPre'][nrnID]
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

def main(meta_data,preEI,preSpikes,savePath,T_MAX=100,activeSoma=False,propofol=0):

    # Initialize neuron morphology
    nrnM = init_neuron(meta_data['mFile'],activeSoma)

    N = T_MAX*16+1
    M = len(meta_data['cellIDs'])
    data = np.zeros([M*5,N])
    # Simulate neurons with given morphology
    for k,nrnID in enumerate(meta_data['cellIDs']):
        cell = LFPy.Cell(nrnM.sectionList,v_init=-65,celsius=37)
        cell.tstop = T_MAX

        # Add synapses to neuron model
        synParams,synTimes = addsyns(nrnID,meta_data,preEI,preSpikes,propofol)
        syn = list()
        for i in range(len(synParams)):
            syn.append(LFPy.Synapse(cell, **synParams[i]))
            syn[-1].set_spike_times(np.array(synTimes[i]))

        # Simulate
        cdm = LFPy.CurrentDipoleMoment(cell=cell)
        cell.simulate(probes=[cdm],rec_somav=True,rec_vmem=False)

        t = cell.tvec.reshape([-1,1])
        v = cell.somav.reshape([-1,1])
        data[5*k:5*(k+1),:] = np.concatenate((t,cdm.data.T,v),axis=1).T
        cell.strip_hoc_objects()

    mType = str.split(meta_data['mFile'],'\\')[-1][:-4]
    saveFile = savePath+'\\'+ mType
    np.save(saveFile, data)


if __name__ == "__main__":
    pars = sys.argv
    n = len(pars)
    postnetwork = pars[1]
    spikingFile = pars[2]
    savePath = pars[3]
    T_MAX = int(pars[4])
    activeSoma = (pars[5]=="true")
    propofol = float(pars[6])

    preEI,preSpikes = load_spike_times(spikingFile)

    conDir = os.path.join(postnetwork,'json')
    conFiles = os.listdir(conDir)

    for meta_data_file in conFiles:
        with open(os.path.join(conDir,meta_data_file), 'r') as f:
            meta_data = json.load(f)
        main(meta_data,preEI,preSpikes,savePath,T_MAX,activeSoma,propofol)