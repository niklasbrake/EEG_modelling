addpath('C:\Users\brake\Documents\MATLAB\fmriView')
% Get surface area of each triangle
[sa,X] = network_simulation.getHeadModel;
x0 = [54,24,35];
[~,i0] = min(vecnorm(X.vertices-x0,2,2));
j0 = find(sum(X.faces==i0,2));
c0 = mean(X.vertices(X.faces(j0(1),:),:));
normal = -sa.cortex75K.normals(i0,:);
[x1,x2] = getOrthBasis(normal);
normal = -[0.8529,0.4924,0.1736];
x1 = [0,0,1];
x2 = [-0.4924,0.8529,0];

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
        if(maxColourChange>0.1)
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
    title('Dipole correlation','fontweight','normal','fontsize',7);
    annotation('line','Position',[0.37 0.82+0.02 0.12 0.10]);
    annotation('line','Position',[0.37 0.715+0.02 0.12 -0.095]);
ax = axes('Position',[0.48,0.575,0.24,0.04]);
    xlim([-20,20]);
    ylim([0,1]);
    text(-10,0.2,'10 mm','FontSize',7,'HorizontalAlignment','center','VerticalAlignment','bottom');
    line([-15,-5],[0,0],'color','k','LineWidth',2);
    axis off;



load('C:\Users\brake\Documents\temp\cortical_area.mat')
N = 16e9;
dMids = 0.5*(rValues(2:end)+rValues(1:end-1));
nrnCount = mean(diff(A),2)*200000;
nrnCount(end) = N-sum(nrnCount(1:end-1));
corr_kernel = @(d) exp(-d.^2/10);
rho_bar = sum(corr_kernel(dMids).*nrnCount)/sum(nrnCount');
SIG_N = @(rho) N+N*(N-1)*rho;

load('E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\model_fits\pre.mat');
load('E:\Research_Projects\004_Propofol\data\simulations\analyzed\_archive\network_criticality_spectra_(s=1).mat')
asynchUnitarySpec = mean(mean(P(:,:,:),3),2);
SIG0 = mean(sum(asynchUnitarySpec*mean(diff(f))));
P0 = interp1(f,asynchUnitarySpec,freq);

red = [0.8000    0.2980    0.0078];
ax = axes('Position',[0.18 0.15 0.23 0.23*8.5/6]);
    plotwitherror(freq,pre,'M','color',[0.6,0.6,0.6],'LineWidth',1);
    plot(freq,P0*50/SIG0,'color',red,'LineWidth',1)
    plot(freq,P0*200/SIG0,'color',red,'LineWidth',1)
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
    txt = text(0.75,0.3,'Model','fontsize',7,'color',red);
    txt = text(1.9,110,'Data','fontsize',7,'color','k');
    gcaformat;
ax = axes('Position',[0.55 0.15 0.23 0.23*8.5/6]);
    dscale = [0,10.^linspace(-2,0,991),1:0.01:10];
    t = [0,10.^linspace(-3,-1,481),0.1:0.001:1];
    t = t(:);
    tt = zeros(length(dscale),length(t));
    for i = 1:length(dscale)
        corr_kernel2 = @(d) exp(-d.^2/dscale(i));
        rho_bar = sum(corr_kernel2(dMids).*nrnCount)/sum(nrnCount');
        tt(i,:) = SIG0*SIG_N(t*rho_bar);
    end
    [XX,YY] = meshgrid(t,dscale);
    surf(XX,YY,tt*0,log10(tt),'LineStyle','none');
    line([0,0.4],[0,0],[0,0],'color','k','LineWidth',0.75)
    line([0,0],[0,10],[0,0],'color','k','LineWidth',0.75)
    view([0,90]);
    hold on;
    contour(XX,YY,log10(tt),log10([50,200]),'color',red,'linewidth',1)
    xlim([0,0.4])
    ylim([0,10]);
    gcaformat
    axis xy
    ylabel('\sigma^2 (mm)')
    xl = xlabel('\rho_{max}');
    xl.Position(2) = -2.13;
    CM = clrsPT.iridescent(1e3);
    CM = interp1(linspace(0,1,1e3),CM,linspace(0,0.98,1e3).^2,'nearest');
    colormap(ax,CM);
    CB = colorbar('location','eastoutside');
    CB.Label.String = 'Total EEG power (uV^2)';
    CB.Label.FontSize=7;
    set(gca,'CLim',[-1,3]);
    CB.Ticks = [-1:3];
    CB.TickLabels = {'0.1','1','10','100','1000'};
    CB.Position = [ax.Position(1)+0.25,ax.Position(2),0.02,ax.Position(4)];
labelpanel(0.07,0.92,'a',true);
labelpanel(0.07,0.48,'b',true);
labelpanel(0.44,0.48,'c',true);

figureNB(5,5);
subplot(1,2,1);
    dscale = [0,10.^linspace(-2,0,991),1.01:0.01:10];
    t = 1;
    tt = zeros(length(dscale),length(t));
    for i = 1:length(dscale)
        corr_kernel2 = @(d) exp(-d.^2/dscale(i));
        rho_bar = sum(corr_kernel2(dMids).*nrnCount)/sum(nrnCount');
        tt(i) = SIG0*SIG_N(t*rho_bar);
    end
    plot(dscale,tt,'LineWidth',1);
    hold on;
    line([1e-1,10],[50,50],'color','r','LineWidth',1);
    line([1e-1,10],[200,200],'color','r','LineWidth',1);
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([1e-1,10]);
    ylim([1,1e4])
    xlabel('\sigma^2 (mm)')
    ylabel('Total EEG power (uV^2)')
    sigma_max = interp1(tt',dscale,200);
    disp(['Since \rho_{max}<1, then \sigma^2 must be greater than ' num2str(sigma_max,2) ' mm.']);
    title(['\rho_{max} = 1']);
    gcaformat;
subplot(1,2,2);
    t = [0,10.^linspace(-3,-1,481),0.1001:0.001:1];
    t = t(:);
    corr_kernel2 = @(d) exp(-d.^2/4);
    rho_bar = sum(corr_kernel2(dMids).*nrnCount)/sum(nrnCount');
    tt = SIG0*SIG_N(t*rho_bar);
    plot(t,tt,'LineWidth',1);
    hold on;
    line([1e-2,1],[50,50],'color','r','LineWidth',1);
    line([1e-2,1],[200,200],'color','r','LineWidth',1);
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([1e-2,1]);
    ylim([1,1e4])
    xlabel('\rho_{max}')
    ylabel('Total EEG power (uV^2)')
    rho_max = interp1(tt,t,200);
    title(['\sigma^2 = 4 mm']);
    disp(['If \sigma^2 = 4 mm, then \rho_{max} should be about ' num2str(rho_max,2) '.']);
    gcaformat;

figureNB(10,5);
subplot(1,2,1);
    plot(rValues,mean(A,2),'.-','MarkerSize',10,'LineWidth',1);
    hold on;
    plot(rValues,mean(B,2),'.-','MarkerSize',10,'LineWidth',1);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlabel('Distance (mm)')
    ylabel('Cumulative area (mm^2)')
    gcaformat
subplot(1,2,2);
    dMids = 0.5*(rValues(2:end)+rValues(1:end-1));
    plot(dMids,mean(diff(A),2),'.-','MarkerSize',10,'LineWidth',1);
    hold on;
    plot(dMids,mean(diff(B),2),'.-','MarkerSize',10,'LineWidth',1);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlabel('Distance (mm)')
    ylabel('Area (mm^2)')
    gcaformat

