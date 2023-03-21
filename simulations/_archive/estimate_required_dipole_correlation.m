
[sa,X] = network_simulation.getHeadModel;

X.vertices = sa.cortex75K.vc;
X.faces= sa.cortex75K.tri;

% Increase grid size of cortex mesh three-fold
Mv = size(X.vertices,1);
M = size(X.faces,1);
V3 = [X.vertices;zeros(M,3)];
F3 = zeros(M*3,3);
p2 = [0,0,0];
I = eye(3);
for i = 1:M
    waitbar(i/M);
    iF = X.faces(i,:);
    V3(Mv+i,:) = mean(X.vertices(iF,:));
    F3(3*(i-1)+(1:3),:) = iF-I.*iF + (Mv+i)*I;
end
X.vertices = V3; 
X.faces = F3; 

% Compute area of each of the triangular faces
M = size(X.faces,1);
A = zeros(M,1);
for i = 1:M
    waitbar(i/M)
    iF = X.faces(i,:);
    p = [X.vertices(iF(1),:); ...
            X.vertices(iF(2),:); ...
            X.vertices(iF(3),:)];
    v1 = p(1,:)-p(2,:);
    v2 = p(1,:)-p(3,:);
    N1 = norm(v1,2);
    N2 = norm(v2,2);
    CosTheta = dot(v1,v2)/(N1*N2);
    A(i) = 0.5*N1*N2*sqrt(1-CosTheta.^2);
end

dThresh = 10.^(linspace(-1,log10(20),60));

% Get each of the three vertices of the mesh triangles faces
F1 = X.vertices(X.faces(:,1),:);
F2 = X.vertices(X.faces(:,2),:);
F3 = X.vertices(X.faces(:,3),:);

% For every point in the cortex, compute the total surface area
% within a radius of dThresh
M = 500;
Ar = zeros(length(dThresh),M);
for kk = 1:M
    waitbar(kk/M);
    i0 = randi(size(X.faces,1));
    p = X.vertices(X.faces(i0,:),:);
    v1 = p(2,:)-p(1,:);
    v2 = p(3,:)-p(1,:);
    u = rand(2,1);
    if(sum(u)>1)
        u = 1-u;
    end
    x0 = p(1,:) + v1*u(1) + v2*u(2);

    % Which of the three vertices are within the radius? Weight the
    % full area of the triangle based on how many vertices are
    % within radius (quick and dirty)
    d1 = vecnorm(F1-x0,2,2);
    d2 = vecnorm(F2-x0,2,2);
    d3 = vecnorm(F3-x0,2,2);
    for i = 1:length(dThresh)
        i1 = find(d2<=dThresh(i));
        idcs1 = find(d1<dThresh(i));
        idcs2 = find(d2<dThresh(i));
        idcs3 = find(d3<dThresh(i));
        Ar(i,kk) = (sum(A(idcs1))+sum(A(idcs2))+sum(A(idcs3)))/3;
    end
end

% Neuron count within a radius is the surface area times
% 100,000 neurons/mmm^2
N = 25e9;
dMids = 0.5*(dThresh(2:end)+dThresh(1:end-1));
nrnCount = mean(diff(Ar),2)*100000;
nrnCount(end) = N-sum(nrnCount(1:end-1));
corr_kernel = @(d) exp(-d.^2/10);

figureNB;
    plot(corr_kernel(dMids'),nrnCount/N,'LineWidth',1);
    set(gca,'yscale','log')
    ylabel('% neurons')
    xlabel('Pairwise correlation')
    gcaformat;

% Average pairwise correlation based on the :
%   - distribution of neuron densities within a certain radius of %     each point in cortex
%   - coupling kernel as a function of radius (normalized, so that
%     the correlation of nearby neurons is 1)
rho_bar = sum(corr_kernel(dMids').*nrnCount)/sum(nrnCount')



load('E:\Research_Projects\004_Propofol\Modelling\neuron_simulations\data\summary_data\asynch_untiary_spectrum.mat');
% asynchUnitarySpec = P;

SIG0 = mean(sum(asynchUnitarySpec*mean(diff(f))));
P0 = interp1(f,asynchUnitarySpec,freq);
load('E:\Research_Projects\004_Propofol\Modelling\data_fitting\data\pre.mat');

SIG_N = @(rho) N+N*(N-1)*rho;

figureNB;
% subplot(1,2,1);
    plotwitherror(freq,0.1*P0*SIG_N(rho_bar),'Q','color','r');
%     set(gca,'xscale','log');
%     set(gca,'yscale','log');
%     xlim([1,300]);
%     xticks([1,10,100]);
%     xticklabels([1,10,100]);
%     xlabel('Frequency (Hz)');
%     ylabel(['PSD (' char(956) 'V^2/Hz)'])
%     ylim([1e-2,1e2])
% subplot(1,2,2);
    plotwitherror(freq,pre,'Q','color',[0.6,0.6,0.6],'LineWidth',1);
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([1,300]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xlabel('Frequency (Hz)');
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    ylim([1e-2,1e2])
gcaformat(gcf);