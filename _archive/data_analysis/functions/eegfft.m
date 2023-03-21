function [freq,time,psd] = eegfft(t,eeg,windowSize,overlap);
%EEGFFT compute spectrogram of EEG signal using multitaper method.
%	[freq,time,psd] = EEGFFT(t,eeg,windowSize,overlap) returns the
%	frequence and time vectors as the first two outputs, and returns
%	the power at each frequency and time point as a matrix (psd).
%	The input t is the time vector for the EEG signal, windowSize is
% 	the length of time over which the EEG signal is segmented and the
% 	power spectrum is computed. The input overlap determines how much
% 	overlap is used between each window.
%
% 	Example:
% 		[freq,time,psd] = eegfft(t,eeg,2,1.9);
% 		imagesc(time,freq,log10(psd));

	Fs = round(1/mean(diff(t))); % sampling frequency
	T = windowSize; % window length in seconds
	iwindow = -floor(Fs*T/2)+1:floor(T/2*Fs);
	L = length(iwindow);
	f = Fs * linspace(0,1/2,L/2+1);
	f = f(2:end);
	minT = floor(Fs*T/2);
	maxT = size(eeg,1)-floor(Fs*T/2);
	difT = ceil(Fs*(T-overlap));
	rachTs = minT:difT:maxT;
	hmw = hamming(L);
	% K = zeros(L/2,length(rachTs));
	K = [];
	thet = K;
	for i = 1:length(rachTs);
		temp = eeg(rachTs(i)+iwindow);
		% if(any(isnan(temp)))
		% 	K(:,i) = nan;
		% 	continue;
		% end
		[K(:,i),f] = pmtm(temp,2,[],Fs); % 3 tapers
		% K(:,i) = 2*temp(2:L/2+1);
	end
	freq = f(2:end);
	time = t(rachTs(:));
	psd = K(2:end,:);