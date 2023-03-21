function [mFiles,synapseFiles,N] = import_postsyanptic_network(networkPath)
    fid = fopen(fullfile(networkPath,'mTypes.txt'),'r');
    temp = textscan(fid,'%s'); 
    mFiles = temp{1};

    synDir = fullfile(networkPath,'connections');
    F = dir(synDir); F = F(3:end);
    n = length(F);
    synapseFiles = cell(n,1);
    for i = 1:n
        synapseFiles{i} = fullfile(synDir,F(i).name);
    end

    N = csvread(fullfile(networkPath,'TotalSynapseCount.txt'));
end