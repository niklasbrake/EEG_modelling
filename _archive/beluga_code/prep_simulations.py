import json
import sys
import os
import numpy as np
import csv
from scipy.stats import poisson
from math import ceil

def import_postsyanptic_network(networkPath):
    file = networkPath + '/mTypes.txt'
    mFiles = list()
    with open(file, newline='') as csvfile:
         syns = csv.reader(csvfile, delimiter='\n')
         for row in syns:
            mFiles.append(row[0])
    file = networkPath + '/connections.csv'
    connections = list()
    with open(file, newline='') as csvfile:
         syns = csv.reader(csvfile, delimiter=',')
         for row in syns:
            connections.append([int(x) for x in row])
    connections = np.array(connections)
    connections[:,1:] = connections[:,1:]-1
    return mFiles,connections

def load_spike_times(spikingFile):
    """ Load spike times from the file specified
    by the variable spikingFile """
    preEI = list()
    preSpikes = list()
    with open(spikingFile, newline='') as csvfile:
         spikeTimes = csv.reader(csvfile, delimiter=' ')
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
         syns = csv.reader(csvfile, delimiter=' ')
         for row in syns:
            temp = row[0].split(',')
            synSeg.append(int(temp[0])-1)
            synPre.append(int(temp[1])-1)
    return synSeg,synPre

def getSynParams(mFile,synapseFiles):
    synSeg = dict()
    synPre = dict()
    IDs = list()
    for k,file in enumerate(synapseFiles):
        nrnID = str.split(file,'/')[-1][:-4]
        synSeg[nrnID],synPre[nrnID] = load_syanpse_locations(file)
        IDs.append(nrnID)
    return synSeg,synPre,IDs


def getSynParams(mFile,synapseFiles):
    synSeg = dict()
    synPre = dict()
    IDs = list()
    for k,file in enumerate(synapseFiles):
        nrnID = str.split(file,'/')[-1][:-4]
        synSeg[nrnID],synPre[nrnID] = load_syanpse_locations(file)
        IDs.append(nrnID)
    return synSeg,synPre,IDs

if __name__ == "__main__":
    pars = sys.argv
    n = len(pars)
    networkPath = pars[1]

    mFiles,connections = import_postsyanptic_network(networkPath)

    mTypes = dict()
    for i,m in enumerate(mFiles):
        i += 1
        if m not in mTypes.keys():
            mTypes[m] = dict({'mFile':m,'cellIDs':list(),'synSeg':dict(),'synPre':dict()})

        mTypes[m]['cellIDs'].append(str(i))
        idcs = np.argwhere(connections[:,0] == i)
        mTypes[m]['synSeg'][i] = connections[idcs,1][:,0].tolist()
        mTypes[m]['synPre'][i] = connections[idcs,2][:,0].tolist()

    if not os.path.exists(networkPath+'/json'):
        os.mkdir(networkPath+'/json')
    for m in mTypes.keys():
        mType = str.split(m,'\\')[-1][:-4]
        with open(networkPath+'/json/meta_data_'+mType+'.json', 'w') as f:
            json.dump(mTypes[m], f)

