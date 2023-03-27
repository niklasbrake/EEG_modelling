load('E:\Research_Projects\005_Aperiodic_EEG\data\simulations\axon_propogation\test_run\long_time_long_space.mat')

% Axoplasmic resisitivty (Ohm um)
Ra = par.node.elec.pas.axres.value*1e6;

% Get section type (node or internode?)
iNode = zeros(size(MEMBRANE_POTENTIAL,2),1);
iNode(nodes) = 1;

% Get comparment lengths
dx_intn = par.intn.seg.geo.length.value.ref;
dx_node = par.node.seg.geo.length.value.ref;
L = iNode*dx_node + (1-iNode)*dx_intn;
SPACE_VECTOR = cumsum(L)*1e-3;

% Get comparment cross-sectional areas
a_intn = pi*par.intn.geo.diam.value.ref.^2/4; % Internodal axon cross sectional area (um)
a_node = pi*par.node.geo.diam.value.ref.^2/4; % Nodal axon cross sectional area (um)
A = iNode*a_node + (1-iNode)*a_intn;


% Current between middle of compartment i to middle of compartment j
x_mids = 0.5*(SPACE_VECTOR(2:end)+SPACE_VECTOR(1:end-1));
R = L(1:end-1)*Ra./(2*A(1:end-1)) + L(2:end)*Ra./(2*A(2:end));
dx = 0.5*(L(1:end-1)+L(2:end));

% Get delta V
dV = -diff(MEMBRANE_POTENTIAL,1,2);

% Compute axial current
Ia = dV./R(:)';
Ia = Ia*1e6; % mA -> nA

% Compute current dipoles
Q = Ia.*dx(:)'; % Dipole (um*nA = fA m)
Q = Q*1e-3; % fA m -> pA m

% Get total axon dipole
apDsingle = sum(Q,2);
