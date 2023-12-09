import sys
import csv
import numpy as np
import scipy as sp
from scipy import optimize
import matplotlib.pyplot as plt
import math

def readdata(filename):
    with open(filename, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=' ', quotechar='|')
        fs = []
        psd = []
        for idx,row in enumerate(reader):
            temp = row[0].split(',')
            fs.append(float(temp[0]))
            psd.append(list(map(float,temp[1:])))
    psd = np.array(psd)
    fs = np.array(fs)
    return fs, psd

def savedata(data,filename):
    with open(filename,'w',newline='') as csvfile:
        wrt = csv.writer(csvfile, delimiter=',')
        if type(data[0]) is list:
            for ps in data:
                if(len(ps)==0):
                    wrt.writerow('0')
                else:
                   wrt.writerow(ps)
        else:
            wrt.writerow(data)

def preparedata(x,y):
    ynew = np.log10(y)
    return np.squeeze(x),np.squeeze(ynew)

def scaleFrequency(f):
    return (2*np.pi*f)**2

def eq1(fScaled, *params):
    tauI, tauE, A1, A2 = params[0]
    ys = np.zeros_like(fScaled)
    x2 = 10**A1*np.reciprocal(1+fScaled*tauI**2) + 10**A2*np.reciprocal(1+fScaled*tauE**2)
    ys = ys + np.log10(x2)
    return ys

def full_model_eq1(xs,*params):
    return eq1(scaleFrequency(xs),params[0][:4]) + spectral_peaks(xs,params[0][4:])


def eq5(fScaled, *params):
    tau1, offset, mag = params[0]
    ys = np.zeros_like(fScaled)
    ys = ys + mag + np.log10(np.exp(offset) + tau1*np.reciprocal(1+fScaled*tau1**2))
    return ys

def full_model_eq5(xs,*params):
    return eq5(scaleFrequency(xs),params[0][:3]) + spectral_peaks(xs,params[0][3:])

def eq6(fScaled, *params):
    tau1,tau2,offset,mag = params[0]
    tau2 = 4e-3;
    ys = np.zeros_like(fScaled)
    x2 = (tau1-tau2)**2 / ((1+tau1**2*fScaled)*(1+tau2**2*fScaled));
    ys = ys + mag + np.log10(np.exp(offset)+x2);
    return ys

def full_model_eq6(xs,*params):
    return eq6(scaleFrequency(xs),params[0][:4]) + spectral_peaks(xs,params[0][4:])

def spectral_peaks(xs, *params):
    ys = np.zeros_like(xs)
    params = params[0]
    for ii in range(0, len(params), 3):
        ctr, hgt, wid = params[ii:ii+3]
        ys = ys + hgt * np.exp(-(xs-ctr)**2 / (2*wid**2))
    return ys

def objective(params, model_func, data):
    xd, yd = data
    ym = model_func(xd,params)
    r = np.mean(np.abs(yd - ym))
    return r

def main(f,p,nPeaks=3,fitType='eq6',sp_ap=list(),sp_p=list()):

    # Periodic parameter bounds, assuming the following peaks
    #   Peak 1: delta rhythm
    #   Peak 2: alpha rhythm
    #   Peak 3: beta rhythm
    lb_p = [0,0,0.2,6,0,0.6,15,0,1,2,0,0.4]
    ub_p = [4,4,3,15,4,4,40,3,10,5,2,2]
    lb_p = lb_p[:3*nPeaks]
    ub_p = ub_p[:3*nPeaks]

    # Aperiodic parameter bounds
    if(fitType=='eq6'):
        lb_ap = [7e-3,3.9e-3,-21,1]
        ub_ap = [75e-3,4.1e-3,-7,5]
        full_model = full_model_eq6
        model_func = eq6
    elif(fitType=='eq1'):
        lb_ap = [4e-3,1e-3,-3,-3]
        ub_ap = [75e-3,20e-3,1,5]
        full_model = full_model_eq1
        model_func = eq1
    elif(fitType=='eq5'):
        lb_ap = [7e-3,-20,1]
        ub_ap = [75e-3,0,5]
        full_model = full_model_eq5
        model_func = eq5

    K = len(lb_ap)

    startpoint = sp_ap
    if(len(sp_ap)==0):
        if(fitType=='eq5'):
            sp_ap = [17e-3,-10.5,4.2]
        elif(fitType=='eq6'):
            sp_ap = [17e-3,4e-3,-10.5,4.2]
        elif(fitType=='eq1'):
            sp_ap = [20e-3,3e-3,0,2]
        sp_p = [0.5,2,1.5,8,0.3,1,22,0.1,4,4,0,1]
        startpoint = sp_ap + sp_p[:3*nPeaks]
    elif(len(sp_ap)==K):
        sp_p = [0.5,2,1.5,8,0.3,1,22,0.1,4,4,0,1]
        startpoint = sp_ap + sp_p[:3*nPeaks]

    # Bounds for full model
    lb = lb_ap + lb_p
    ub = ub_ap + ub_p

    # Get initial parameter values by fitting each component seperately
    [x_data,y_data] = preparedata(f,p)

    if len(y_data.shape)==1:
        # Fit just aperiodic component
        results1 = sp.optimize.least_squares(objective,startpoint,bounds=(lb,ub),args = [full_model,[x_data, y_data]])
        parsSave = results1.x
        if(nPeaks>0):
            # Fit peaks, using results1 as initial conditions for aperiodic component
            yDentrended = y_data-model_func(scaleFrequency(x_data),results1.x[:K])
            results2 = sp.optimize.least_squares(objective,results1.x[K:],bounds=(lb_p,ub_p),args = [spectral_peaks,[x_data, yDentrended]])

            pars0 = np.concatenate((results1.x[:K],results2.x),0)
            results = sp.optimize.least_squares(objective,pars0,bounds=(lb,ub),args = [full_model,[x_data, y_data]])
            parsSave = results.x
    else:
        # If input is a matrix, use fit to previous timepoint as initial condition
        results1 = sp.optimize.least_squares(objective,startpoint,bounds=(lb,ub),args = [full_model,[x_data, y_data[:,0]]])
        pars0 = results1.x

        m = y_data.shape[1]
        parsSave = np.zeros([m,len(pars0)])
        for i in range(m):
            for j,par in enumerate(pars0):
                if(par>ub[j]):
                    pars0[j] = ub[j]
                if(par<lb[j]):
                    pars0[j] = lb[j]
            results = sp.optimize.least_squares(objective,pars0,bounds=(lb,ub),args = [full_model,[x_data, y_data[:,i]]])
            parsSave[i,:] = results.x
    return parsSave

if __name__ == "__main__":
    filename = sys.argv[1]
    fitType = sys.argv[2]
    nPeaks = int(sys.argv[3])
    startPointFile = sys.argv[4]
    f,p = readdata(filename)

    if(len(startPointFile)>0):
        with open(startPointFile, newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            startpoint = list()
            for row in reader:
                for value in row:
                    startpoint.append(float(value))
    else:
        startpoint = list()

    pars = main(f,p,nPeaks,fitType,startpoint)
    saveName = filename[:-4] + '_params.csv'
    savedata(pars.tolist(),saveName)
    print(saveName,end='')

