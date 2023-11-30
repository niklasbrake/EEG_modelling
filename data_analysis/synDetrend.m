function [params,synFun,full_model] = synDetrend(f,P,nPeaks,fitType,startPoint)
    myPath = strrep(fileparts(mfilename('fullpath')),'\','/');
    if(nargin<3)
        nPeaks = 3;
    elseif(~isnumeric(nPeaks) || nPeaks>4)
        error('nPeaks must be an integer less than 3')
    end
    if(nargin<4)
        fitType = 'exp2';
    elseif(~strcmp(fitType,'exp2') && ~strcmp(fitType,'lorenz') && ~strcmp(fitType,'unilorenz') && ~strcmp(fitType,'avalanches'))
        error('fitType must be either ''exp2'' or ''lorenz'' or ''unilorenz'' or ''avalanches''');
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
