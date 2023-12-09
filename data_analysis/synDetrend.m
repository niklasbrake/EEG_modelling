function [params,synFun,full_model] = synDetrend(f,P,nPeaks,fitType,startPoint)
% MATLAB wrapper for detrending.py
    myPath = strrep(fileparts(mfilename('fullpath')),'\','/');
    if(nargin<3)
        nPeaks = 3;
    elseif(~isnumeric(nPeaks) || nPeaks>4)
        error('nPeaks must be an integer less than 3')
    end
    if(nargin<4)
        fitType = 'eq6';
    elseif(~strcmp(fitType,'eq6') && ~strcmp(fitType,'eq1') && ~strcmp(fitType,'eq5') && ~strcmp(fitType,'avalanches'))
        error('fitType must be either ''eq6'' or ''eq1'' or ''eq5'' or ''avalanches''');
    end
    if(nargin<5)
        startPointsFile = '""';
    else
        startPointsFile = strrep(fullfile(myPath,'startPoint.csv'),'\','/');
        csvwrite(startPointsFile,startPoint);
    end
    fileName = 'temp';
    file0 = strrep(fullfile(myPath,[fileName '.csv']),'\','/');
    pyFunction = fullfile(myPath,'detrending.py');
    f = f(:);
    if(size(P,1)==1)
        P = P(:);
    end
    csvwrite(file0,[f,P]);
    cmd = ['python ' pyFunction ' ' file0 ' ' fitType ' ' int2str(nPeaks) ' ' startPointsFile];
    [err,file] = system(cmd);
    if(err)
        error(file)
    end
    params = csvread(file);
    delete(file0);
    if(nargin==5)
        delete(startPointsFile);
    end
    [full_model,synFun] = fittingmodel(fitType);
end
