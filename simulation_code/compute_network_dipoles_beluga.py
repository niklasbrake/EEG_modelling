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

def addsyns(nrnID,connection_data,preEI,preSpikes,parameters):
    # SYANPSE_FAILURE_RATE = 0.7
    SYANPSE_FAILURE_RATE = 0

    synSeg = connection_data['synSeg'][nrnID]
    synPre = connection_data['synPre'][nrnID]
    synTimes = list()
    synParams = list()
    # Get excitatory syanpse locations
    for i in range(len(synSeg)):
        ts = preSpikes[synPre[i]]
        spikeSelection = np.random.binomial(len(ts),1-SYANPSE_FAILURE_RATE)
        ts = np.random.choice(ts,spikeSelection,replace=False).tolist()
        if(preEI[synPre[i]]):
            params = parameters["iSynParams"].copy()
        else:
            params = parameters["eSynParams"].copy()
        params['idx'] = synSeg[i]
        if(len(ts)>0):
            synParams.append(params)
            synTimes.append(ts)

    return synParams,synTimes

def main(parameters,connection_data,preEI,preSpikes,savePath,T_MAX=100):

    # Initialize neuron morphology
    nrnM = init_neuron(connection_data['mFile'],parameters)

    N = T_MAX*16+1
    M = len(connection_data['cellIDs'])
    data = np.zeros([M*5,N])
    # Simulate neurons with given morphology
    for k,nrnID in enumerate(connection_data['cellIDs']):
        cell = LFPy.Cell(nrnM.sectionList,v_init=parameters["biophys_pars"]['pas_mem_pars']['erev_leak'],celsius=37)
        cell.tstop = T_MAX

        # Add synapses to neuron model
        synParams,synTimes = addsyns(nrnID,connection_data,preEI,preSpikes,parameters)
        syn = list()
        for i in range(len(synParams)):
            syn.append(LFPy.Synapse(cell, **synParams[i]))
            syn[-1].set_spike_times(np.array(synTimes[i]))

        # Simulate
        cdm = LFPy.CurrentDipoleMoment(cell=cell)
        cell.simulate(probes=[cdm],rec_somav=True,rec_vmem=True)

        t = cell.tvec.reshape([-1,1])
        v = cell.somav.reshape([-1,1])
        data[5*k:5*(k+1),:] = np.concatenate((t,cdm.data.T,v),axis=1).T
        # multi_dipoles, dipole_locs = cell.get_multi_current_dipole_moments()
        cell.strip_hoc_objects()


    mType = str.split(connection_data['mFile'],'/')[-1][:-4]
    saveFile = savePath+'/'+ mType
    np.save(saveFile, data)

    # from scipy.io import savemat
    # fm = saveFile + '_multi_dipoles.mat'
    # savemat(fm, {'dipoles':multi_dipoles,'x':dipole_locs,'time':t})

if __name__ == "__main__":
    pars = sys.argv
    n = len(pars)
    postnetwork = pars[1]
    spikingFile = pars[2]
    savePath = pars[3]
    T_MAX = int(pars[4])
    if(len(pars)<6):
        parFile = '_parameters.json'
    else:
        parFile = pars[5]

    with open(os.path.join(savePath,parFile), 'r') as f:
        parameters = json.load(f)

    preEI,preSpikes = load_spike_times(spikingFile)

    conDir = os.path.join(postnetwork,'json')
    conFiles = os.listdir(conDir)

    for connection_data_file in conFiles:
        with open(os.path.join(conDir,connection_data_file), 'r') as f:
            connection_data = json.load(f)
        main(parameters,connection_data,preEI,preSpikes,savePath,T_MAX)