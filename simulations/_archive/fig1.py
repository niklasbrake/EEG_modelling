import numpy as np
from neuron import h
from matplotlib import pyplot as plt
import LFPy

def morph2Segs(cellID):
    data = np.loadtxt(cellID)
    temp = np.array(data[:,3])
    data[:,3] = data[:,4]
    data[:,4] = temp
    arbours = getChildDendrite(1,data)
    segs = list()
    cons = list()
    for i in arbours:
        segs,cons = getSegments(np.array([i]),1,-1,segs,cons,data)

    endPoints = [[x[i] for i in (0,-1)] for x in segs]
    c0 = np.array([data[i[0]-1,2:6] for i in endPoints])
    c1 = np.array([data[i[1]-1,2:6] for i in endPoints])
    X = np.concatenate((c0,c1),1)
    # plotMorphology(segs,data)
    return X,cons,segs,data

def getSegments(node,nodepar,segpar,segs,cons,data):
    i = len(segs)
    segs.append(np.append(nodepar,advanceSegment(node,data)))
    # segs.append(np.array(advanceSegment(node,data)))
    cons.append([segpar,i])
    child = getChildDendrite(segs[-1][-1],data).tolist()
    if(len(child)>0):
        segs,cons = getSegments([child[0]],segs[i][-1],i,segs,cons,data)
        for j in range(1,len(child)):
            segs,cons = getSegments([child[j]],segs[i][-1],i,segs,cons,data)
    return segs,cons

def advanceSegment(branch,data):
    node = branch[-1]
    childIdcs = getChildDendrite(node,data)
    if(len(childIdcs)==1):
        if(len(branch)==1):
            branch = np.append(branch,childIdcs)
            branch = advanceSegment(branch,data)
        else:
            L = np.sum(np.linalg.norm(np.diff(data[branch-1,2:5],1,0),2,1)) # DEFINE MAX SEGMENT LENGTH TO BE 20 UM
            if(L<=20):
                branch = np.append(branch,childIdcs)
                branch = advanceSegment(branch,data)
    elif(len(branch)>1):
        branch = branch[:-1]
    return branch

def getChildDendrite(node,data):
    isChild = data[node-1:,-1]==node
    isDendrite = (data[node-1:,1]==3)+(data[node-1:,1]==4)
    return node+np.argwhere(isDendrite*isChild).flatten()

def plotMorphology(segs,data):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    for i in range(len(segs)):
        X = data[segs[i]-1,2:5]
        plt.plot(X[:,0],X[:,1],X[:,2])
    axisEqual3D(ax)
    plt.show()

def axisEqual3D(ax):
    extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    sz = extents[:,1] - extents[:,0]
    centers = np.mean(extents, axis=1)
    maxsize = max(abs(sz))
    r = maxsize/2
    for ctr, dim in zip(centers, 'xyz'):
        getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)
