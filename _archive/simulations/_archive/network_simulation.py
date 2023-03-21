import numpy as np

def simulatespikes(tmax,M,m,rateOnly=False):
    if(m>1):
        error('Branching index must be less than 1.')
    tmax = 5
    m = 0.9
    M = 1000
    tmax = tmax+2
    dt = 4e-3
    t = np.arange(0,tmax,dt)
    N = len(t)

    ME = int(np.floor(0.85*M))
    MI = int(M-ME)

    k=4
    h = 1*dt*ME*(1-m)

    # Generate mean firing rate using critical branching process
    B = np.zeros([N,1])
    exN = np.random.poisson(h,[N,1])
    for i in range(N-1):
        if(B[i]):
            count = np.random.binomial(B[i][0]*k,m/k)
        else:
            count = 0
        B[i+1] = count+exN[i]

    iRemove = np.argwhere(t<0)
    B = np.delete(B,iRemove)
    t = np.delete(t,iRemove)-2

    if(rateOnly):
        return t,B

    numspikes = int(np.sum(B))
    spikeTime = np.zeros(5*numspikes)
    neuronIDs = np.zeros(5*numspikes)
    lRatio = 1
    j = 0
    for i in range(len(t)):
        nEx = np.random.poisson(B[i])
        spikeTime[j:j+nEx] = t[i] + np.random.random(nEx)*dt - dt/2
        neuronIDs[j:j+nEx] = np.random.choice(ME,nEx,replace=False)
        j+=nEx

        nIn = np.random.poisson(lRatio*B[i])
        spikeTime[j:j+nIn] = t[i] + np.random.random(nIn)*dt - dt/2
        neuronIDs[j:j+nIn] = ME+np.random.choice(MI,nIn,replace=False)
        j+=nIn
    spikeTime = spikeTime[:j]
    neuronIDs = neuronIDs[:j]
    ei = np.concatenate((np.zeros(ME),np.ones(MI)),0)
    I = np.argsort(neuronIDs)
    neuronIDs = neuronIDs[I]
    spikeTime = 1e3*spikeTime[I]

    return t,B,neuronIDs,spikeTime,ei