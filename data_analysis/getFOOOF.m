function [OOF,fun,peaks,rsquared] = getFOOOF(f,P,withKnee)

fileName = 'temp';
myPath = strrep(fileparts(mfilename('fullpath')),'\','/');
file0 = strrep(fullfile(myPath,[fileName '.csv']),'\','/');
csvwrite(file0,[f,P]);

% Analyze with FOOOF in python. Outputs to .csv files
if(withKnee)
	cmd = ['python ' fullfile(myPath,'analyzePSDpeaks.py') ' ' file0 ' ' int2str(1)];
	fun = @(x,p) 10^p(1) * (p(2) + x.^p(3)).^(-1);
else
	cmd = ['python ' fullfile(myPath,'analyzePSDpeaks.py') ' ' file0 ' ' int2str(2)];
	fun = @(x,p) 10^p(1)./x.^p(2);
end
[err,prints] = system(cmd);
if(err)
	error(prints)
end


file1 = fullfile(myPath,[fileName '_peakFreq.csv']);
file2 = fullfile(myPath,[fileName '_peakWidth.csv']);
file3 = fullfile(myPath,[fileName '_peakHeight.csv']);
file4 = fullfile(myPath,[fileName '_OOF.csv']);
file5 = fullfile(myPath,[fileName '_rsquared.csv']);


peakFreq = csvread(file1);
peakWidth =  csvread(file2);
peakHeight = csvread(file3);
OOF = csvread(file4);
rsquared = csvread(file5);

delete(file0);
delete(file1);
delete(file2);
delete(file3);
delete(file4);
delete(file5);
peaks.freq = peakFreq;
peaks.width = peakWidth;
peaks.Height = peakHeight;