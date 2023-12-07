dataFolder = 'E:\Research_Projects\004_Propofol\manuscript\Version3\Data';
load(fullfile(dataFolder,'simulations_synapse_dipole_orientation.mat'));

dMag = squeeze(max(vecnorm(dipoles,2,2)));
dMag./max(dMag)*3;

dipole_moment = squeeze(nanmedian(dipoles./vecnorm(dipoles,2,2)))';
dipole_moment = [dipole_moment(:,1),dipole_moment(:,3),-dipole_moment(:,2)];
dipole_moment(find(E_or_I==0),:) = -dipole_moment(find(E_or_I==0),:);

clrs = [230, 25, 75;
    60, 180, 75];

[theta0,phi0] = cart2sph(synapse_position(:,1),synapse_position(:,2),synapse_position(:,3));
[theta1,phi1] = cart2sph(dipole_moment(:,1),dipole_moment(:,2),dipole_moment(:,3));
synapse_position_pol = [theta0,phi0];
PD_pol = [theta1,phi1];

% gIdcs = cellfun(@(x)isempty(strfind(x,'E')),{network.neurons.mType});
gIdcs = 1-E_or_I;
clrs = [0,0,1;1,0,0];
clrs = clrs(gIdcs+1,:);

% Remove points close to the corners of the plot (to avoid variation crossing periodic boundary conditions being displayed)
idcs1 = find(~and(PD_pol(:,1)<-0.85*pi,synapse_position_pol(:,1)>0.85*pi));
idcs2 = find(~and(PD_pol(:,2)<-0.85*pi/2,synapse_position_pol(:,2)>0.85*pi/2));
idcs = intersect(idcs1,idcs2);
% Room for label
idcs1 = find(~and(PD_pol(:,1)>0.85*pi,synapse_position_pol(:,1)<0));
idcs2 = find(~and(PD_pol(:,2)>0.75*pi/2,synapse_position_pol(:,2)<0.5));
idcs = intersect(idcs,intersect(idcs1,idcs2));
% Make data sparse for display purposes
idcs = idcs(randperm(size(idcs,1),400));


fig = figureNB(12.9,3.5);
ax(1) = axes('Position',[0.08, 0.30, 0.15, 0.55])
ax(2) = axes('Position',[0.29, 0.30, 0.15, 0.55])
for i = 1:2
    axes(ax(i));
    scatter(synapse_position_pol(idcs,i),PD_pol(idcs,i),dMag(idcs),'k','filled'); hold on;
    FT = fitlm(synapse_position_pol(idcs,i),PD_pol(idcs,i),'Weights',dMag(idcs))
end
axes(ax(1));
    gcaformat;
    title('Azimuth','FontSize',7,'FontWeight','normal');
    ylabel('Dipole moment')
    xlim([-pi,pi]); xticks([-pi,0,pi]);
    ylim([-pi,pi]); yticks([-pi,0,pi]);
    xticklabels({sprintf('%c%c',char(45),char(960)),'0',sprintf('%c',char(960))})
    yticklabels({sprintf('%c%c',char(45),char(960)),'0',sprintf('%c',char(960))})
    line([-pi,pi],[-pi,pi],'color','r');
axes(ax(2));
    gcaformat;
    title('Elevation','FontSize',7,'FontWeight','normal');
    xl = xlabel('Synapse location (rel. soma)');
    xl.Position = [-2.2,-2.5,-1];
    xlim([-pi,pi]/2); xticks([-pi/2,0,pi/2]);
    ylim([-pi,pi]/2); yticks([-pi/2,0,pi/2]);
    xticklabels({sprintf('%c%c/2',char(45),char(960)),'0',sprintf('%c/2',char(960))})
    yticklabels({sprintf('%c%c/2',char(45),char(960)),'0',sprintf('%c/2',char(960))})
    line([-pi,pi]/2,[-pi,pi]/2,'color','r');


load('E:\Research_Projects\004_Propofol\data\simulations\raw\synthetic_spikes\synthetic_spikes_analyzed.mat');
shuffle = load('E:\Research_Projects\004_Propofol\data\simulations\raw\synthetic_spikes\synthetic_spikes_space_analyzed.mat');

for i = 1:size(dipoles,4)
    C(i,1) = corr(dipoles(:,1,1,i),dipoles(:,1,2,i));
    C(i,2) = corr(dipoles(:,2,1,i),dipoles(:,2,2,i));
    C(i,3) = corr(dipoles(:,3,1,i),dipoles(:,3,2,i));
end

for i = 1:size(dipoles,4)
    C_spaceshuffle(i,1) = corr(shuffle.dipoles(:,1,1,i),shuffle.dipoles(:,1,2,i));
    C_spaceshuffle(i,2) = corr(shuffle.dipoles(:,2,1,i),shuffle.dipoles(:,2,2,i));
    C_spaceshuffle(i,3) = corr(shuffle.dipoles(:,3,1,i),shuffle.dipoles(:,3,2,i));
end


cr = reshape(nanmean(C,2),[6,11])';
cr_shuffle = reshape(nanmean(C_spaceshuffle,2),[6,11])';
rRange = 10.^linspace(-4,-1,6);

m0 = nanmean(cr);
s0 = stderror(cr);

m1 = nanmean(cr_shuffle);
s1 = stderror(cr_shuffle);

axes('Position',[0.55, 0.30, 0.15, 0.55])
    plot(rRange,m0,'.-k','LineWidth',1,'MarkerSize',10)
    hold on;
    line([rRange;rRange],[m0-1.96*s0;m0+1.96*s0],'color','k','LineWidth',1)
    plot(rRange,m1,'.-','color',[0.6,0.6,0.6],'LineWidth',1,'MarkerSize',10)
    hold on;
    line([rRange;rRange],[m1-1.96*s1;m1+1.96*s1],'color',[0.6,0.6,0.6],'LineWidth',1)
    set(gca,'xscale','log');
    xlim([1e-4,0.0252])
    ylabel('Dipole correlation')
    xlabel('Synapse correlation (R_{max})')
    gcaformat


load('E:\Research_Projects\004_Propofol\data\simulations\raw\synthetic_spikes\synthetic_spikes_analyzed.mat');
idcs = find(R==2);
dipoles = dipoles(:,:,:,idcs);
dipoles = reshape(dipoles,160001,3,[]);
dp = resample(dipoles,1e3,16e3);

[sa,X] = network_simulation_beluga.getHeadModel;
location = randi(size(sa.cortex75K.vc,1),1e3,1);
eegPSD = zeros(8193,1);
for k = 1:length(location)
    waitbar(k/length(location));
    eeg = network_simulation_beluga.getEEG(dp,sa,location(k));
    [psd,f] = pmtm(detrend(eeg,'constant'),2,[],1e3);
    eegPSD = eegPSD + mean(psd,2);
end
eegPSD = eegPSD/length(location);

% Get baseline EEG spectrum from propofol cohort
dataFolder = 'E:\Research_Projects\004_Propofol\manuscript\Nature Communications\_final_submission\_data'
data = load(fullfile(dataFolder,'electrode2_Cz.mat'));
data.baseline = squeeze(nanmedian(data.psd(:,data.tRescaled<-1,:),2));

axes('Position',[0.81 0.30, 0.15, 0.55])
    plotwitherror(data.freq,data.baseline,'M','color',[0.5,0.5,0.5]);
    plot(f,eegPSD*1e15,'color','r');
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.5,100]);
    xlabel('Frequency (Hz)');
    xticks([1,10,100])
    xticklabels([1,10,100])
    ylabel('PSD (uV^2/Hz)');
    xax = get(gca,'xaxis');
    xax.TickLabelRotation=0;
    gcaformat;



