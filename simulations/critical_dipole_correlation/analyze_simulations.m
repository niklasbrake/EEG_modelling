function initialize_network(masterFolder,S)
    addpath('/lustre04/scratch/nbrake/code/simulation_code');

    count = 1;
    C = zeros(100,100,3);
    for arrayID = 1:10
        for j = 1:10
            i = 10*(arrayID-1)+j;
            folder = fullfile(masterFolder,['run_' num2str(i,'%.3d')],'simulation');
            for rep = 1:100
                file = fullfile(folder,['simulation_S=' sprintf('%0.1f',S) '_rep' int2str(rep) '.mat']);
                if(~exist(file))
                    C(i,rep,:) = nan;
                    continue;
                end
                load(fullfile(file));
                C(i,rep,1) = corr(dipoles(:,1,1),dipoles(:,1,2));
                C(i,rep,2) = corr(dipoles(:,2,1),dipoles(:,2,2));
                C(i,rep,3) = corr(dipoles(:,3,1),dipoles(:,3,2));
            end
        end
    end
    time = time(1:16:end);
    % save(fullfile(masterFolder,'dipole_correlation_(s=0).mat'),'C');
    save(fullfile(masterFolder,['dipole_correlation_S=' sprintf('%0.1f',S) '.mat']),'C');
end
