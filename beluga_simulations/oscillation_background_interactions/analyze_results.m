% addpath('/lustre04/scratch/nbrake/code/simulation_code');

folder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\peak_interactions';
% folder = '/lustre04/scratch/nbrake/data/simulations/peak_interactions';
% importData(folder);

load(fullfile(folder,'analyzed_data.mat'))
fid = fopen('C:\Users\brake\Documents\GitHub\EEG_modelling\simulation_code\mod_files\default_parameters.json');
str = char(fread(fid,inf)');
fclose(fid);
defaultPars = jsondecode(str);
[f,P,par] = organizeSimulationData(psd,EI,mValue,tau,defaultPars);


function [f,P,par] = organizeSimulationData(psd,EI,mValue,tau,defaultPars)
    defaultEI = defaultPars.eCellParams.firingRate/defaultPars.iCellParams.firingRate;
    defaultEI = 1e-3*round(1e3*defaultEI);

    % Segment data by parameters
    % Tau changes
    idcsTau = find(tau~=defaultPars.iSynParams.tau2);
    P_tau = cellfun(@(x)mean(x,2),psd(idcsTau),'UniformOutput',false);
    P_tau = cat(2,P_tau{:});
    tauAx = tau(idcsTau);

    % EI changes
    idcsEI = find(EI~=defaultEI);
    P_EI = cellfun(@(x)mean(x,2),psd(idcsEI),'UniformOutput',false);
    P_EI = cat(2,P_EI{:});
    EIAx = EI(idcsEI);

    % Branching changes
    idcsM = find(mValue~=0);
    P_m = cellfun(@(x)mean(x,2),psd(idcsM),'UniformOutput',false);
    P_m = cat(2,P_m{:});
    mAx = mValue(idcsM);

    % Default parameter sets
    idcsDefault = setdiff(1:length(psd),union(union(idcsTau,idcsEI),idcsM));
    P_default = mean(cat(2,psd{idcsDefault}),2);

    % Supplement each changing vector with the default
    P_tau = cat(2,P_tau,P_default);
    tauAx = [tauAx,defaultPars.iSynParams.tau2];
    [tauAx,I] = sort(tauAx);
    P_tau = P_tau(:,I);

    P_EI = cat(2,P_EI,P_default);
    EIAx = [EIAx,defaultEI];
    [EIAx,I] = sort(EIAx,'descend');
    P_EI = P_EI(:,I);

    P_m = cat(2,P_m,P_default);
    mAx = [mAx,0];
    [mAx,I] = sort(mAx);
    P_m = P_m(:,I);

    P.tau = P_tau;
    P.EI = P_EI;
    P.m = P_m;
    par.tau = tauAx;
    par.m = mAx;
    par.EI = EIAx;

    f = 0.1:0.1:1e3;
end


function importData(folder)
    [sa,X] = network_simulation_beluga.getHeadModel;
    F = dir(folder);
    goodIdcs = find(and(~strcmp({F(:).name},'..'),~strcmp({F(:).name},'.')));
    F = F(goodIdcs);
    goodIdcs = find([F(:).isdir]);
    F = F(goodIdcs);

    locations = randi(size(sa.cortex75K.vc,1),20,1); % Random location

    for j = 1:length(F)
        folder2 = fullfile(folder,F(j).name);
        F2 = dir(folder2);
        F2 = F2(3:end);
        dp = [];
        for i = 1:length(F2)
            try
                folder3 = fullfile(folder2,F2(i).name,'simulation');
                data = load(fullfile(folder3,'simulation_data.mat'));
                dipoles = resample(data.dipoles(2:end,:),2e3,16e3);
                dp(:,:,i) = dipoles;
            catch
            end
        end
        psd{j} = zeros(10e3,size(dp,3));
        for k = 1:length(locations)
            eeg = network_simulation_beluga.getEEG(dp,sa,locations(k));
            psd{j} = psd{j} + mypmtm(detrend(eeg,'constant'),2e3,10)/length(locations);
        end
    end

    % Set up parameter values
    for j = 1:length(F)
        folder2 = fullfile(folder,F(j).name);
        F2 = dir(folder2);
        F2 = F2(3:end);
        dp = [];
        i=1;
        folder3 = fullfile(folder2,F2(i).name,'simulation');
        load(fullfile(folder2,F2(i).name,'model.mat'));
        mValue(j) = network.branchNo;
        fid = fopen(fullfile(folder3,'_parameters.json'));
        str = char(fread(fid,inf)');
        fclose(fid);
        pars = jsondecode(str);
        EI(j) = pars.eCellParams.firingRate/pars.iCellParams.firingRate;
        tau(j) = pars.iSynParams.tau2;
    end
    mValue = mValue;
    EI = 1e-3*round(1e3*EI);
    tau = tau;

    save(fullfile(folder,'analyzed_data.mat'),'psd','EI','tau','mValue');
end