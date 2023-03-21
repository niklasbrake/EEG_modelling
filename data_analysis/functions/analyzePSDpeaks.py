import csv
import sys
import numpy as np
from fooof import FOOOF

def main(filename,withKnee):
    # Read .csv file
    with open(filename, newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=' ', quotechar='|')
        fs = []
        psd = []
        for idx,row in enumerate(spamreader):
            temp = row[0].split(',')
            fs.append(float(temp[0]))
            psd.append(list(map(float,temp[1:])))
    # Turn into numpy arrays
    psd = np.array(psd)
    fs = np.array(fs)
    # Run fooof on each time average from this patient
    if(withKnee  == "1"):
        fm = FOOOF(peak_width_limits=[2.0, 15.0], max_n_peaks=3, aperiodic_mode='knee')
    else:
        fm = FOOOF(peak_width_limits=[2.0, 15.0], max_n_peaks=3, aperiodic_mode='fixed')
    peakFreq = []
    peakHeight = []
    OOF = []
    peakWidth = []
    rsquare = []
    for i in range(0,psd.shape[1]):
        fm.fit(fs,psd[:,i])
        peakFreq.append(fm.get_results().peak_params[:,0])
        peakWidth.append(fm.get_results().peak_params[:,1])
        peakHeight.append(fm.get_results().peak_params[:,2])
        OOF.append(fm.aperiodic_params_)
        rsquare.append(fm.r_squared_)
    # Save peak position and height into .csv files
    with open(filename[:-4] + '_peakFreq.csv','w',newline='') as csvfile:
        wrt = csv.writer(csvfile, delimiter=',')
        for peaks in peakFreq:
            if(len(peaks)==0):
                wrt.writerow('0')
            else:
                wrt.writerow(peaks)
    with open(filename[:-4] + '_peakHeight.csv','w',newline='') as csvfile:
        wrt = csv.writer(csvfile, delimiter=',')
        for peaks in peakHeight:
            if(len(peaks)==0):
                wrt.writerow('0')
            else:
                wrt.writerow(peaks)
    with open(filename[:-4] + '_OOF.csv','w',newline='') as csvfile:
        wrt = csv.writer(csvfile, delimiter=',')
        for ps in OOF:
            if(len(ps)==0):
                wrt.writerow('0')
            else:
               wrt.writerow(ps)
    with open(filename[:-4] + '_peakWidth.csv','w',newline='') as csvfile:
        wrt = csv.writer(csvfile, delimiter=',')
        for ps in peakWidth:
            if(len(ps)==0):
                wrt.writerow('0')
            else:
                wrt.writerow(ps)
    with open(filename[:-4] + '_rsquared.csv','w',newline='') as csvfile:
        wrt = csv.writer(csvfile, delimiter=',')
        wrt.writerow(rsquare)

if __name__ == "__main__":
    # Take input arguments from matlab
    filename = sys.argv[1]
    withKnee = sys.argv[2]
    main(filename,withKnee)
