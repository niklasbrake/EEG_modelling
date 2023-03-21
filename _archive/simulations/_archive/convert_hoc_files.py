from hoc2swc import hoc2swc
from os import listdir
folder = 'E:/Research_Projects/004_Propofol/Modelling/neuron_simulations/data/cortical_column_Hagen'

for f in listdir(folder +'/hoc'):
    hoc2swc(folder + '/hoc/' + f,folder + '/swc/' + f[:-4] + '.swc')

filepath = 'E:/Research_Projects/004_Propofol/Modelling/neuron_simulations/data/cortical_column_Hagen/swc/'

for neuronID in os.listdir(filepath):
    file = filepath + neuronID
    pts3d,connections,segs,morphData = getMorphoSegments.morph2Segs(file)
    badSegs = list()
    for i in range(len(segs)):
        for j in range(len(segs[i])-1):
            idx1 = segs[i][j]
            idx2 = segs[i][j+1]
            L = np.sum((morphData[idx1-1,2:5]-morphData[idx2-1,2:5])**2)
            if(L<1e-8):
                badSegs.append([idx1,idx2])
    for i in range(len(badSegs)):
        idcs = np.argwhere(morphData[:,0]==badSegs[i][1])
        morphData[idcs,3]+=0.1 
        # for j in idcs:
        #     morphData[j,-1] = badSegs[0][0]
        # idcs = np.argwhere(morphData[:,0]==badSegs[i][1])
        # morphData=np.delete(morphData,idcs[0],0)
    np.savetxt(filepath+neuronID, morphData, fmt=' '.join(['%i'] + ['%i'] + ['%1.4f'] + ['%1.4f'] + ['%1.4f'] + ['%1.4f'] + ['%i']))