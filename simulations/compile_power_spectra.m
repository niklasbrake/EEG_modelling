folder = '/lustre04/scratch/nbrake/data/simulations/raw/parameter_search';
F = dir(folder); F = F(3:end);

P = [];
dgam = [];
dlam = [];
dtau = [];
tau = [];
for i = 1:length(F)
    filename = fullfile(folder,F(i).name,'power_spectra.mat');
    if(exist(filename)~=0);
        data = load(filename);
        m = size(data.P,2);
        if(m<200)
            data.P = cat(2,data.P,nan(size(data.P,1),200-m));
        end
        P = cat(3,P,data.P);
        dgam = cat(2,dgam,data.dgam);
        dlam = cat(2,dlam,data.dlam);
        dtau = cat(2,dtau,data.dtau);
        tau = cat(2,tau,data.tau);
    end
end
f = data.f;
save(fullfile(folder,'power_spectra_all.mat'),'P','dgam','dlam','dtau','tau','f');