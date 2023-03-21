import numpy as np
from neuron import h
import LFPy

def morph2Segs(cellID):
    data = np.loadtxt(cellID)
    arbours = list()
    somaPoints = np.argwhere(data[:,1]==1)
    for i in somaPoints:
        temp = get_child(int(data[i,0]),data).tolist()
        if(len(temp)>0):
            arbours.append((i,temp))
    segs = list()
    cons = list()
    for i in range(len(arbours)):
        for j in arbours[i][1]:
            segs,cons = get_segments(np.array([j]),arbours[i][0],-1,segs,cons,data)

    endPoints = [[x[i] for i in (0,-1)] for x in segs]
    c0 = np.array([data[i[0]-1,2:6] for i in endPoints])
    c1 = np.array([data[i[1]-1,2:6] for i in endPoints])
    compType0 = np.array([data[i[0]-1,1,None] for i in endPoints])
    compType1 = np.array([data[i[1]-1,1,None] for i in endPoints])
    pts3d = np.concatenate((c0,c1,compType0,compType1),1)
    # plot_morphology(segs,data)
    return pts3d,cons,segs,data

def get_segments(node,parentNode,parentSegment,segs,cons,data):
    i = len(segs)
    segs.append(np.append(parentNode,advance_segment(node,data)))
    # segs.append(np.array(advance_segment(node,data)))
    cons.append([parentSegment,i])
    child = get_child(segs[-1][-1],data).tolist()
    if(len(child)>0):
        segs,cons = get_segments([child[0]],segs[i][-1],i,segs,cons,data)
        for j in range(1,len(child)):
            segs,cons = get_segments([child[j]],segs[i][-1],i,segs,cons,data)
    return segs,cons

def advance_segment(branch,data):
    node = branch[-1]
    childIdcs = get_child(node,data)
    if(len(childIdcs)==1):
        if(len(branch)==1):
            branch = np.append(branch,childIdcs)
            branch = advance_segment(branch,data)
        else:
            L = np.sum(np.linalg.norm(np.diff(data[branch-1,2:5],1,0),2,1)) # DEFINE MAX SEGMENT LENGTH TO BE 20 UM
            if(L<=20):
                branch = np.append(branch,childIdcs)
                branch = advance_segment(branch,data)
    elif(len(branch)>1):
        branch = branch[:-1]
    return branch

def get_child(node,data):
    isChild = data[node-1:,-1]==node
    isDendrite = (data[node-1:,1]==3)+(data[node-1:,1]==4)
    return node+np.argwhere(isDendrite*isChild).flatten()

def axis_equal(ax):
    extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    sz = extents[:,1] - extents[:,0]
    centers = np.mean(extents, axis=1)
    maxsize = max(abs(sz))
    r = maxsize/2
    for ctr, dim in zip(centers, 'xyz'):
        getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)
