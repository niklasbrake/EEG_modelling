function initialize_network(masterFolder)
    addpath('/lustre04/scratch/nbrake/code/simulation_code');

    count = 1;
    C = zeros(100,100,3);
    dp = zeros(2001,3,10000);
    for arrayID = 1:10
        for j = 1:10
            i = 10*(arrayID-1)+j;
            folder = fullfile(masterFolder,['run_' num2str(i,'%.3d')]);
            for rep = 1:100
                % saveFolder = fullfile(folder,['simulation' int2str(rep)]);
                saveFolder = fullfile(folder,['simulation' int2str(rep)]);
                file = fullfile(saveFolder,'simulation_data.mat');
                if(exist(file))
                    load(file);
                    dp(:,:,100*(i-1)+rep) = resample(sum(dipoles,3),1e3,16e3);
                else
                    dp(:,:,100*(i-1)+rep) = nan;
                end
            end
        end
    end
    time = time(1:16:end);
    dipoles = dp;
    save(fullfile(masterFolder,'dipoles.mat'),'time','dipoles');
end
