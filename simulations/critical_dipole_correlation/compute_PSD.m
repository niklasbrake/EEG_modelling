function initialize_network(masterFolder)

addpath('/lustre04/scratch/nbrake/code/simulation_code');

[sa,X] = network_simulation_beluga.getHeadModel;
locations = sa.cortex2K.in_from_cortex75K;

% dp = zeros(2001,3,100);
% for arrayID = 1:10
%     for j = 1:10
%         i = 10*(arrayID-1)+j;
%         folder = fullfile(masterFolder,['run_' num2str(i,'%.3d')]);
%         for rep = 1:100
%             % saveFolder = fullfile(folder,['simulation' int2str(rep)]);
%             saveFolder = fullfile(folder,['simulation' int2str(rep)]);
%             load(fullfile(saveFolder,'simulation_data.mat'));
%             dp(:,:,rep) = dp(:,:,rep) + resample(sum(dipoles,3),1e3,16e3);
%         end
%     end
% end
% time = time(1:16:end);
% dipoles = dp;
load(fullfile(masterFolder,'dipoles.mat'),'dipoles','time');

m=10;
for i = 1:10
    dp = dipoles(:,:,m*(i-1)+1:m*i);
    eeg = zeros(size(dp,1),m*length(locations));
    G = zeros(1,m*length(locations));
    for k = 1:length(locations)
        eeg(:,m*(k-1)+1:m*k) = network_simulation_beluga.getEEG(dp,sa,locations(k));
        G(m*(k-1)+1:m*k) = 1:m;
    end
    psd = mypmtm(detrend(eeg,'constant'),1e3,2);
    Ptemp = splitapply(@(x) mean(x,2),psd,G);
    if(i==1)
        P = nan(size(psd,1),size(dp,3));
    end
    P(:,m*(i-1)+1:m*i) = Ptemp;
end

f = [0.5:0.5:500];
save(fullfile(masterFolder,'PSD.mat'),'f','P');

function psd = mypmtm(xin,fs,bins_per_hz)
% Taken from https://www.mathworks.com/matlabcentral/answers/267759-how-to-adapt-pmtm-for-a-gpuarray
    [m, n] = size(xin);%m=length of signal, n=# of signals
    k=fs*bins_per_hz;
    nfft = 2^nextpow2(m+k-1);%- Length for power-of-two fft.
    [E,V] = dpss(m,4);
    g=gpuDevice;
    s=(m+nfft)*32*length(V);%how many bytes will be needed for ea. signal
    ne=floor(g.AvailableMemory/s/2);%number of signals that can be processed at once with available memory
    indx=[0:ne:n,n];%number of iterations that will be necessary
    psd=zeros(k/2,n);%initialize output
      for i=1:length(indx)-1
        x=gpuArray(xin(:,1+indx(i):indx(i+1)));
        w = exp(-1i .* 2 .* pi ./ k);
        x=x.*permute(E,[1 3 2]); %apply dpss windows
        %------- Premultiply data.
        kk = ( (-m+1):max(k-1,m-1) )';
        kk = (kk .^ 2) ./ 2;
        ww = w .^ (kk);   % <----- Chirp filter is 1./ww
        nn = (0:(m-1))';
        x = x .* ww(m+nn);
        %------- Fast convolution via FFT.
        x = fft(  x, nfft );
        fv = fft( 1 ./ ww(1:(k-1+m)), nfft );   % <----- Chirp filter.
        x = x .* fv;
        x  = ifft( x );
        %------- Final multiply.
        x = x( m:(m+k-1), : , :) .* ww( m:(m+k-1) );
        x = abs(x).^2;
        x=x.*permute(V,[2 3 1])/length(V);%'eigen' method of estimating average
        x=sum(x,3);
        x=x(2:end/2+1,:)./fs;
        psd(:,1+indx(i):indx(i+1))=gather(x);
      end
  end