function sampled_parameters = sample_parameter_space(N)
    % Define distributions of parameters
    lambda_E = @(p) icdf('logn',p,log(0.5),1);
    lambda_I = @(p) icdf('logn',p,log(2.5),1);
    tauE = @(p) icdf('unif',p,1,3.5);
    tauI = @(p) icdf('unif',p,5,20);
    erev = @(p) icdf('unif',p,-75,-45);
    gleak = @(p) 10.^icdf('unif',p,-5,log10(0.005));

    % Distribution of morphologies
    neuronTypes = {'L23E_oi24rpy1';'L23I_oi38lbc1';'L23I_oi38lbc1';'L4E_53rpy1';'L4E_j7_L4stellate';'L4E_j7_L4stellate';'L4I_oi26rbc1';'L4I_oi26rbc1';'L5E_oi15rpy4';'L5E_j4a';'L5I_oi15rbc1';'L5I_oi15rbc1';'L6E_51_2a_CNG';'L6E_oi15rpy4';'L6I_oi15rbc1';'L6I_oi15rbc1'};
    nrnAbundance = [26.8,3.2,4.3,9.5,9.5,9.5,5.6,1.5,4.9,1.3,0.6,0.8,14,4.6,1.9,1.9];
    mTypeCount = accumarray(findgroups(neuronTypes),nrnAbundance);
    mTypeP = cumsum(mTypeCount)/sum(mTypeCount);
    mType = @(p) interp1(mTypeP,1:11,p,'next','extrap');

    samplePars = @(P) [lambda_E(P(:,1)), lambda_I(P(:,2)), tauE(P(:,3)), ...
                        tauI(P(:,4)), erev(P(:,5)), gleak(P(:,6)), mType(P(:,7))];

    LH = lhsdesign(100,7);
    PARS = samplePars(LH);

    sampled_parameters = array2table(PARS);
    sampled_parameters.Properties.VariableNames = {'eCellParams.firingRate','iCellParams.firingRate','eSynParams.tau2','iSynParams.tau2','biophys_pars.pas_mem_pars.erev_leak','biophys_pars.pas_mem_pars.g_leak','mType'};