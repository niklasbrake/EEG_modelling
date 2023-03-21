import sys
import os
sys.path.append('C:/Users/brake/Documents/pyNB')
sys.path.append('C:/Users/brake/Documents/GitHub/aperiodic_EEG_modelling/simulations/functions')
sys.path.append('C:/Users/brake/Documents/GitHub/aperiodic_EEG_modelling/simulations/_archive')
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

class cortical_neuron():
    def __init__(self,mFile,activeSoma=False):
        AXIAL_RESISTANCE        = 100
        EREV_LEAK               = -63       # Masoli et al
        G_LEAK                  = 0.0011

        pts3d,connections,segs,morphData = getMorphoSegments.morph2Segs(mFile)
        load_mechanisms('C:/Users/brake/Documents/GitHub/aperiodic_EEG_modelling/simulations/_archive/_mod')
        # List to record the type of neuron segment:
        # soma = 0, basal dendrite = 3, apical dendrite = 4
        secType = dict()
        # Generate NEURON list of sections
        self.sectionList = h.SectionList()
        mType = str.split(mFile,'\\')[-1][:-4]
        cnabar, gkbar, shift, shift2 = self.conductanceLU(mType)
        # Add soma to NEURON model
        self.soma = h.Section(name='soma')
        secType['soma'] = 1
        self.soma.insert('Leak')
        self.soma.e_Leak = EREV_LEAK
        self.soma.gmax_Leak = G_LEAK
        if activeSoma:
            self.soma.insert('inaT')
            self.soma.shift_inaT = shift
            self.soma.vtraub_inaT = -68+shift2
            self.soma.gnabar_inaT = cnabar
            self.soma.ena = 60
            self.soma.insert('ikdT')
            self.soma.gkbar_ikdT = gkbar
            self.soma.ek = -80
            self.soma.insert('imZ')
            self.soma.gkbar_imZ = 5*1e-5
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

    def conductanceLU(self,mType):
        shift2 = 0
        if mType == 'L23E_oi24rpy1':
            shift = -15+1
            gnabar = 0.6*1.8
            gkbar = 0.5*1.8*0.6
        elif mType == 'L23I_oi38lbc1':
            shift = -15
            gnabar = 0.6*1.5
            gkbar = 0.5*1.5
        elif mType == 'L4E_53rpy1':
            shift = -15
            shift2 = 15
            gnabar = 0.6
            gkbar = 0.5
        elif mType == 'L4E_j7_L4stellate':
            shift = -15
            gnabar = 0.6*1.5
            gkbar = 0.5*1.5
        elif mType == 'L4I_oi26rbc1':
            shift = -15
            gnabar = 0.6/1.6
            gkbar = 0.5/1.6
        elif mType == 'L5E_j4a':
            shift = -10+3
            gnabar = 0.6/2.5*2
            gkbar = 0.5/2.5*2
        elif mType == 'L5E_oi15rpy4':
            shift = -7
            gnabar = 0.6*5/10*1.5
            gkbar = 0.5*5/10*1.5
        elif mType == 'L5I_oi15rbc1':
            shift = -10
            gnabar = 0.6*1.1
            gkbar = 0.5*1.1
        elif mType == 'L6E_51_2a_CNG':
            shift = -7
            gnabar = 0.6/1.5*0.6
            gkbar = 0.5/1.5*0.6
        elif mType == 'L6E_oi15rpy4':
            shift = -15
            gnabar = 0.6*1.8
            gkbar = 0.5*1.8
        elif mType == 'L6I_oi15rbc1':
            shift = -5
            gnabar = 0.6
            gkbar = 0.5
        else:
            raise ValueError('mType did not match value in lookup table')
        return gnabar, gkbar, shift, shift2