load('/home/nbrake/data/resources/head_models/sa_nyhead.mat');
idcs = sa.cortex2K.in_from_cortex75K;

masterPath = 'E:\Research_Projects\004_Propofol\data\simulations\raw\synapse_embedding_2';
F = dir(masterPath);
F = F(3:end);

m = sort([0.999,0.997,0.99,0.98,0.95,0.86,0.63,0]);

C = zeros(length(F),1);
P = zeros(2e3,length(F2),length(F));
for i = 1:length(F)
    folder = fullfile(masterPath,F(i).name);
    F2 = dir(folder);
    F2 = F2(3:end);
    for j = 1:length(F2)
        for rep = 1:10
            folder2 = fullfile(folder,F2(j).name,['LFPy_' int2str(rep)]);
            load(fullfile(folder2,'simulation_data.mat'));
            M = size(dipoles,3);
            for k = 1:M
                temp = resample(dipoles(2:end-0.2*16e-3,:,k),2e3,16e3);
                if(k==1)
                    dp = zeros(size(temp,1),3,M);
                end
                dp(:,:,k) = temp;
            end
            q = gpuArray(dp);

            % idcs = sa.cortex2K.in_from_cortex75K;
            [f,temp] = expectedEEGspectrum(dp,1,sa,idcs);
            P(:,j,i) = P(:,j,i)+temp;
        end
        P(:,j,i) = P(:,j,i)/10;
    end
end


% clrs = clrsPT.sequential(length(F)+4);
% clrs = clrs(5:end,:);
% figureNB;
% for i = 1:length(F)
%     plotwitherror(f,P(:,:,i),'SE','color',clrs(i,:));
%     hold on;
% end
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% xlim([0.5,150])
% xlabel('Frequency (Hz)');
% ylabel('Power (uV^2/Hz)');
% colormap(clrs);
% CB = colorbar;
% CB.Ticks = linspace(0.5/length(m),1-0.5/length(m),length(m))
% CB.TickLabels = m;
% gcaformat;

function [f,P] = expectedEEGspectrum(dipoles,neuronIDs,sa,idcs)
    q = gpuArray(dipoles(:,:,neuronIDs));
    for j = 1:length(idcs)
        eeg = getEEG(q,sa,idcs(j));
        temp = mypmtm(eeg,2e3,2);
       if(j==1)
            P = zeros(size(temp));
        end
        P = P+temp;
    end
    P = P/length(idcs)*pi/2;
    f = 0.5:0.5:1e3;
end
function psd = mypmtm(xin,fs,bins_per_hz)
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
function W = getEEG(Q,sa,idx)
    % Lead field for Cz recording
    L0 = squeeze(sa.cortex75K.V_fem(49,idx,:))'; % czIDX = 49

    % Orient dipole "up" direction to be normal to cortical surface
    vz = sa.cortex75K.normals(idx,:);
    [vx,vy] = getOrthBasis(vz);
    q = (Q(:,1,:).*vx+Q(:,2,:).*vy+Q(:,3,:).*vz)*1e-12; % nAum -> mAm

    W = squeeze(sum(L0.*q,2)*1e6); % V -> uV
end
function [x1,x2] = getOrthBasis(x0)
    if(prod(size(x0))==3)
        if(x0>0.9)
            x1 = [0,1,0];
        else
            x1 = [1,0,0];
        end
        x1 = x1-x0*sum(x1.*x0);
        x1 = x1/norm(x1,2);
        x2 = cross(x1,x0);
    else
        idcs = x0(:,1)>0.9;
        x1 = zeros(size(x0));
        x1(idcs,2) = 1;
        x1(~idcs,1) = 1;

        x1 = x1-x0.*sum(x1.*x0,2);
        x1 = x1./vecnorm(x1,2,2);
        x2 = cross(x1,x0);
    end
end