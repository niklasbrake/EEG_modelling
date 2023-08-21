from scipy.io import savemat
import numpy as np
import glob
import sys

def main(path,filename,toDelete):
    npzFiles = glob.glob(path + "/*.npy")
    for i,f in enumerate(npzFiles):
        d = np.load(f)
        if(i==0):
            Qx = d[range(1,d.shape[0],5),:]
            Qy = d[range(2,d.shape[0],5),:]
            Qz = d[range(3,d.shape[0],5),:]
            Vm = d[range(4,d.shape[0],5),:]
            t = d[0,:]
        else:
            Qx = np.concatenate((Qx,(d[range(1,d.shape[0],5),:])),axis=0)
            Qy = np.concatenate((Qy,(d[range(2,d.shape[0],5),:])),axis=0)
            Qz = np.concatenate((Qz,(d[range(3,d.shape[0],5),:])),axis=0)
            Vm = np.concatenate((Vm,(d[range(4,d.shape[0],5),:])),axis=0)
            t = d[0,:]

    Qx = np.transpose(np.expand_dims(Qx,2),(1,2,0))
    Qy = np.transpose(np.expand_dims(Qy,2),(1,2,0))
    Qz = np.transpose(np.expand_dims(Qz,2),(1,2,0))
    dipoles = np.concatenate((Qx,Qy,Qz),axis=1)
    fm = path+'/'+filename
    savemat(fm, {'dipoles':dipoles,'time':t.T,'V':Vm.T})
    print('Done (',path, ')')

    if(toDelete):
        import os
        for f in npzFiles:
            os.remove(f)

if __name__ == "__main__":
    pars = sys.argv
    path = pars[1]
    if(len(pars)>2):
        filename = pars[2]
    else:
        filename = 'simulation_data.mat'
    if(len(pars)>3):
        toDelete = (pars[3]=='True')
    else:
        toDelete = False
    main(path,filename,toDelete)

