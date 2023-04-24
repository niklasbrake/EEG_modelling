import sys
import os
import numpy as np
import LFPy
import getMorphoSegments
import csv
from neuron import h, load_mechanisms
from scipy.stats import poisson
from math import ceil
import json


class cortical_neuron():
    mechanismFolder = os.path.join(os.path.dirname(__file__),'mod_files')
    def __init__(self,mFile,parameters):
        pars = parameters['biophys_pars']

        AXIAL_RESISTANCE        = pars['pas_mem_pars']['Ra']
        EREV_LEAK               = pars['pas_mem_pars']['erev_leak']
        G_LEAK                  = pars['pas_mem_pars']['g_leak']

        pts3d,connections,segs,morphData = getMorphoSegments.morph2Segs(mFile)
        load_mechanisms(cortical_neuron.mechanismFolder)

        # List to record the type of neuron segment:
        # soma = 0, basal dendrite = 3, apical dendrite = 4
        secType = dict()
        # Generate NEURON list of sections
        self.sectionList = h.SectionList()
        mType = str.split(mFile,'/')[-1][:-4]
        # Add soma to NEURON model
        self.soma = h.Section(name='soma')
        secType['soma'] = 1
        self.soma.insert('Leak')
        self.soma.e_Leak = EREV_LEAK
        self.soma.gmax_Leak = G_LEAK
        if pars['activeSoma']:
            self.soma.insert('inaT')
            self.soma.shift_inaT = pars['shift'][mType]
            self.soma.vtraub_inaT = -68
            self.soma.gnabar_inaT = pars['gnabar'][mType]
            self.soma.ena = 60
            self.soma.insert('ikdT')
            self.soma.gkbar_ikdT = pars['gkbar'][mType]
            self.soma.ek = -80
            # self.soma.insert('imZ')
            # self.soma.gkbar_imZ = 5*1e-5
        self.soma.push()
        h.pt3dclear()
        i = 0
        xSoma = morphData[morphData[:,1]==1,2:6]
        for x in xSoma:
            h.pt3dadd(x[0], x[1], x[2], x[3])

        self.sectionList.append()
        h.pop_section()
        # Add dendrite segments to NEURON model
        self.dend = list()
        for i in range(len(segs)):
            secName = 'dend'+str(i)
            self.dend.append(h.Section(name=secName))
            self.dend[-1].push()
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
            self.dend[-1].insert('Leak')
            self.dend[-1].e_Leak = EREV_LEAK
            self.dend[-1].gmax_Leak = G_LEAK
            self.dend[-1].Ra = AXIAL_RESISTANCE
            self.sectionList.append()
            h.pop_section()
        # Connect all compartments together
        for i in range(len(connections)):
            i1 = connections[i][0]
            i0 = connections[i][1]
            if(i1==-1):
                self.dend[i0].connect(self.soma,1,0)
            else:
                self.dend[i0].connect(self.dend[i1](1))