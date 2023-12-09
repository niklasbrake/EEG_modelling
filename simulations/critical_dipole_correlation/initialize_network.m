function initialize_network(masterFolder,branchNo,arrayID)
    addpath('/lustre04/scratch/nbrake/code/simulation_code');

    fid = fopen('dyad_morphologies.txt');
    txt = textscan(fid,'%s%s','Delimiter',',');
    mTypes = cat(2,txt{:,1},txt{:,2});

    for j = 1:10
        i = 10*(arrayID-1)+j;
        m_ID(1) = mName2ID(mTypes{i,1});
        m_ID(2) = mName2ID(mTypes{i,2});
        folder = fullfile(masterFolder,['run_' num2str(i,'%.3d')]);
        if(~exist(fullfile(folder,'model.mat')))
            main(folder,m_ID,branchNo)
        end
    end

end

function main(folder,mTypes,branchNo)

    % Initialize network
    network = network_simulation_beluga(folder);

    % Initialize post network
    network = network.initialize_postsynaptic_network(2,mTypes);

    network.spikingFile = fullfile(network.preNetwork,'spikeTimesLong.csv');
    % Presyanptic network parameters
    nPreNeurons = 30e3;
    network.tmax = 40e3;
    network.branchNo = branchNo;
    network.simulatespikes_critplane(nPreNeurons,network.tmax);
    network.save();
end

function m_ID = mName2ID(m_str)
    mTypes = {'L23E_oi24rpy1','L23I_oi38lbc1','L4E_53rpy1','L4E_j7_L4stellate','L4I_oi26rbc1','L5E_j4a','L5E_oi15rpy4','L5I_oi15rbc1','L6E_51_2a_CNG','L6E_oi15rpy4','L6I_oi15rbc1'};
    m_ID = find(strcmp(m_str,mTypes));
end