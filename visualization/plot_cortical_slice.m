% Get surface area of each triangle
x0 = [53.7,23,36.5];
[~,i0] = min(vecnorm(X.vertices-x0,2,2));
j0 = find(sum(X.faces==i0,2));
c0 = mean(X.vertices(X.faces(j0(1),:),:));
normal = -sa.cortex75K.normals(i0,:);
[x1,x2] = getOrthBasis(normal);

% Change coordinate system
V = zeros(size(X.vertices));
V(:,1) = sum((X.vertices-c0).*x2,2);
V(:,2) = sum((X.vertices-c0).*normal,2);
V(:,3) = sum((X.vertices-c0).*x1,2);
idcs = zeros(size(X.faces,1),1);
% Remove faces far away
for i = 1:size(X.faces);
    tri = V(X.faces(i,:),:);
    if(all(vecnorm(tri,Inf,2)<30))
        idcs(i) = 1;
    end
end
F = X.faces(find(idcs),:);


sigma = 4;
pairs = nchoosek(1:3,2);
V0 = V;
F0 = F;
V2 = V;
F2 = [];
newSub = true;
k = 1;
while(newSub)
    newSub = false;
    disp(['round ' int2str(k)])
    for i = 1:size(F0,1)
        tri = V0(F0(i,:),:);
        v_dist = exp(-vecnorm(tri,2,2).^2/sigma);
        maxColourChange = max(abs(diff(v_dist(pairs),1,2)));
        if(maxColourChange>0.05)
            [v1,f1] = subdivide_tri(V2,F0(i,:));
            V2 = [V2;v1];
            newSub = true;
        else
            f1 = F0(i,:);
        end
        F2 = [F2;f1];
    end
    V0 = V2;
    F0 = F2;
    F2 = [];
    k = k+1;
end

X2.vertices = V;
X2.faces = F;
X3.vertices = V0;
X3.faces = F0;

for i = 1:size(X3.faces,1)
    v = X3.vertices(X3.faces(i,:),:);
    v = v-v(1,:);
    v = v(2:end,:)./vecnorm(v(2:end,:),2,2);
    X3.normals(i,:) = -cross(v(1,:),v(2,:));
end

for i = 1:length(X3.vertices)
    fs = find(sum(X3.faces==i,2));
    if(~isempty(fs))
        X3.vec_normals(i,:) = mean(X3.normals(fs,:));
    else
        X3.vec_normals(i,:) = zeros(1,3);
    end
end

% [H,kg] = compute_mesh_curvature(X2.vertices,X2.faces);

sigma = 10;
fig = figure('color','w','units','centimeters');
fig.Position(3:4) = [8.5,6];
ax = axes('Position',[0.11,0.515,0.4,0.5]);
    plot_mesh_brain(X);
    d = vecnorm(X.vertices-c0,2,2);
    C = exp(-d.^2/sigma);
    paint_mesh(C./max(C));
    view([120,10]);
    fix_lighting;
    colormap(ax,clrsPT.sequential(100))
    set(gca,'CLim',[0,1]);

    VW = get(gca,'view');
    [a1,a2,a3] = sph2cart((VW(1)-90)*pi/180,VW(2)/180*pi,1);
    e0 = [a1,a2,a3]';
    e2 = camup';
    e1 = cross(e0,-e2);
    R = 40;
    H = 15;
    x1 = c0(:) + sqrt(2)*R/2*(e1*cos(-3*pi/4)+e2*sin(-3*pi/4))+H*e0;
    x2 = x1(:) + R*(e1*cos(pi/2)+e2*sin(pi/2));
    x3 = x2(:) + R*(e1*cos(0)+e2*sin(0));
    x4 = x3(:) + R*(e1*cos(-pi/2)+e2*sin(-pi/2));
    x = [x1,x2,x3,x4];
    P = patch(x(1,:),x(2,:),x(3,:),'b','LineWidth',1,'FaceColor','none');

ax = axes('Position',[0.455,0.64,0.3,0.3]);
    plot_mesh_brain(X3);
    d = vecnorm(X3.vertices,2,2);
    C = exp(-d.^2/sigma);
    paint_mesh(C/max(C));
    colormap(ax,clrsPT.sequential(100))
    set(gca,'CLim',[0,1]);
    yl = get(gca,'ylim');
    xlim([-20,20]);
    zlim([-20,20]);
    x = [-20,-20,20,20];
    y = [1,1,1,1]*yl(1);
    z = [-20,20,20,-20];
    P = patch(x,y,z,'b','LineWidth',1,'FaceColor','none');
    CB = colorbar;
    CB.Position = [0.74,ax.Position(2),0.02,ax.Position(4)];
    CB.Ticks = [0,1];
    CB.TickLabels = {'0','\rho_{max}'};
    annotation('line','Position',[0.37 0.82+0.02 0.12 0.10]);
    annotation('line','Position',[0.37 0.715+0.02 0.12 -0.095]);
ax = axes('Position',[0.48,0.57,0.24,0.04]);
    xlim([-20,20]);
    ylim([0,1]);
    text(-10,0.2,'10 mm','FontSize',7,'HorizontalAlignment','center','VerticalAlignment','bottom');
    line([-15,-5],[0,0],'color','k','LineWidth',2);
    axis off;

dThresh = 10.^(linspace(-1,log10(20),60));
load('E:\Research_Projects\004_Propofol\data\resources\head_models\surface_area_calculations.mat')
% Neuron count within a radius is the surface area times
% 100,000 neurons/mmm^2
N = 25e9;
dMids = 0.5*(dThresh(2:end)+dThresh(1:end-1));
nrnCount = mean(diff(Ar),2)*100000;
nrnCount(end) = N-sum(nrnCount(1:end-1));
corr_kernel = @(d) exp(-d.^2/10);

% Average pairwise correlation based on the :
%   - distribution of neuron densities within a certain radius of %     each point in cortex
%   - coupling kernel as a function of radius (normalized, so that
%     the correlation of nearby neurons is 1)
corr_kernel = @(d) exp(-d.^2/10);
rho_bar = sum(corr_kernel(dMids').*nrnCount)/sum(nrnCount')


load('E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\model_fits\pre.mat');

load('E:\Research_Projects\004_Propofol\data\simulations\analyzed\network_criticality_spectra_(s=1).mat')

asynchUnitarySpec = mean(mean(P(:,:,:),3),2);
SIG0 = mean(sum(asynchUnitarySpec*mean(diff(f))));
P0 = interp1(f,asynchUnitarySpec,freq);

% asynchUnitarySpec = mean(mean(P(:,1:3,:),3),2);
% SIG0 = mean(sum(asynchUnitarySpec*mean(diff(f))));
% P0 = interp1(f,asynchUnitarySpec,freq);

SIG_N = @(rho) N+N*(N-1)*rho;

ax = axes('Position',[0.18 0.15 0.23 0.23*8.5/6]);
    plotwitherror(freq,pre,'M','color',[0.6,0.6,0.6],'LineWidth',1);
    plot(freq,P0*50/SIG0,'r','LineWidth',1)
    plot(freq,P0*200/SIG0,'r','LineWidth',1)
    % plot(freq,(P0*100/SIG0+P0*100/SIG0)/2,'b','LineWidth',1)
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.5,50]);
    xticks([0.5,5,50])
    xticklabels([0.5,5,50])
    xlabel('Frequency (Hz)');
    ylim([1e-2,2e2])
    yticks([1e-2,1,1e2])
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    gcaformat;
ax = axes('Position',[0.55 0.15 0.23 0.23*8.5/6]);
    dscale = 10.^linspace(-1,1,1e3);
    t = 10.^linspace(-2,log10(01),1e3);
    tt = [];
    for i = 1:length(dscale)
        corr_kernel2 = @(d) exp(-d.^2/dscale(i));
        rho_bar = sum(corr_kernel2(dMids').*nrnCount)/sum(nrnCount');
        tt(:,i) = SIG0*SIG_N(t*rho_bar);
    end
    imagesc((t),(dscale),log10(tt)); hold on;
    contour((t),(dscale),log10(tt),log10([50,200]),'color','r')
    xlim([0,0.4])
    ylim([0,10]);
    gcaformat
    axis xy
    ylabel('\sigma (mm)')
    xl = xlabel('\rho_{max}');
    xl.Position(2) = -2.13;
    colormap(ax,clrsPT.iridescent(1e3))
    CB = colorbar('location','eastoutside');
    CB.Label.String = 'Total EEG power (uV^2)';
    set(gca,'CLim',[-1,3]);
    CB.Ticks = [-1:3];
    CB.TickLabels = {'0.1','1','10','100','1000'};
    CB.Position = [ax.Position(1)+0.25,ax.Position(2),0.02,ax.Position(4)];

labelpanel(0.07,0.92,'a',true);
labelpanel(0.07,0.48,'b',true);
labelpanel(0.44,0.48,'c',true);