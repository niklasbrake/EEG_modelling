from neuron import h
import numpy as np

class cell2comp:
    def __init__(self):
        self.soma = h.Section(name='soma')
        self.soma.nseg = 1
        self.soma.diam = 29.8
        self.soma.cm = 0.77
        self.soma.L =   29.8
        self.soma.Ra = 122

        self.soma.insert('pas')


        self.dend = h.Section(name='dend')
        self.dend.insert('pas')
        self.dend.connect(self.soma,1,0)

        self.syn = h.Exp2Syn(self.soma(0.5))
        self.syn.e = -80

        self.s = h.NetStim()
        self.s.interval = 100
        self.s.start = 20
        self.s.noise = 1

        self.ns = h.NetCon(self.s,self.syn)
        self.ns.weight[0] = 1

        self.t = h.Vector()
        self.t.record(h._ref_t)
        
        self.vm_soma = h.Vector()
        self.vm_soma.record(self.soma(0.5)._ref_v)

        self.vm_dend = h.Vector()
        self.vm_dend.record(self.dend(0.5)._ref_v)