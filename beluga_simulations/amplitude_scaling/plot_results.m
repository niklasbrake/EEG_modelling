psd = [];
% Get surface area of each triangle
[sa,X] = network_simulation_beluga.getHeadModel;
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
    set(ax,'CLim',[0,1]);
    colormap(ax,clrsPT.sequential(100))
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

% Compute expected EEG variance of 16 billion passive neurons
load(fullfile(dataFolder,'anatomy_cortical_pairwise_distance_distribution.mat'));
signed_area = A;
total_area = B;
N = 16e9;
dMids = 0.5*(rValues(2:end)+rValues(1:end-1));
nrnCount = mean(diff(signed_area),2)*200000;
nrnCount(end) = N-sum(nrnCount(1:end-1));
corr_kernel = @(d) exp(-d.^2/10);
rho_bar = sum(corr_kernel(dMids).*nrnCount)/sum(nrnCount');
SIG_N = @(rho) N+N*(N-1)*rho;

% Get baseline EEG spectrum from propofol cohort
data = load(fullfile(dataFolder,'electrode2_Cz.mat'));
data.baseline = squeeze(nanmedian(data.psd(:,data.tRescaled<-1,:),2));

% Get simulated passive spectrum
model = load(fullfile(dataFolder,'simulated_spectra_parameter_sensitivity.mat'));
asynchUnitarySpec = mean(model.psd(:,:,:),2);
SIG0 = sum(asynchUnitarySpec*mean(diff(model.f)));
P0 = interp1(model.f,asynchUnitarySpec,data.freq);

A1 = 30;
A2 = 200;

% Plot
red = [0.8000    0.2980    0.0078];
red = [0,0,0];
ax = axes('Position',[0.18 0.15 0.23 0.23*8.5/6]);
    plotwitherror(data.freq,data.baseline,'M','LineWidth',1,'color',[0.6,0.6,0.6]);
    plot(data.freq,P0*A1/SIG0,'color',red,'LineWidth',1)
    plot(data.freq,P0*A2/SIG0,'color',red,'LineWidth',1)
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
    dscale = [0,10.^linspace(-2,0,991),1:0.01:13];
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
    C = contour(XX,YY,log10(tt),log10([A1,A2]),'color','k','linewidth',1);
    C1 = C(:,2:C(2,1)+1);
    C2 = C(:,C(2,1)+3:end);
    % P = patch([C1(1,:),fliplr(C2(1,:))],[C1(2,:),fliplr(C2(2,:))],'b','LineStyle','none');
    % H = hatchfill(P,'single',45,3,'none');
    % H.Color = red;
    xlim([0,0.2])
    ylim([0,13]);
    gcaformat
    axis xy
    ylabel('\sigma^2 (mm)')
    xl = xlabel('\rho_{max}');
    xl.Position(2) = -2.13;
    CM = clrsPT.iridescent(1e3);
    CM = interp1(linspace(0,1,1e3),CM,linspace(0,0.98,1e3).^2,'nearest');
    % colormap(ax,CM);
    CB = colorbar('location','eastoutside');
    CB.Label.String = 'Total EEG power (uV^2)';
    CB.Label.FontSize=7;
    % set(gca,'CLim',[-1,3]);
    set(gca,'CLim',[0,3])
    CB.Ticks = [-1:3];
    CB.TickLabels = {'0.1','1','10','100','1000'};
    CB.Position = [ax.Position(1)+0.25,ax.Position(2),0.02,ax.Position(4)];
    CM = clrsPT.diverging(1e3);
    c = linspace(0,3,1e3);
    t = linspace(-1,1,1e3);
    cScale = (c-log10(125)).*(c<log10(A1))./(log10(125))+(c-log10(125)).*(c>log10(A2))./(3-log10(125));
    % cScale = (c-log10(A1))./(0.05-(c-log10(A1))).*(c<log10(A1)) + ...
     % (c-log10(A2))./(0.05+(c-log10(200))).*(c>log10(200));

    CM2 = interp1(t,CM,cScale);
    colormap(ax,CM2)


labelpanel(0.07,0.92,'a',true);
labelpanel(0.07,0.48,'b',true);
labelpanel(0.44,0.48,'c',true);





dscale = [0,10.^linspace(-2,0,991),1.01:0.01:13.1];
t = [0,10.^linspace(-3,-1,481),0.101:0.001:1];
t = t(:);
tt = zeros(length(dscale),length(t));
for i = 1:length(dscale)
    corr_kernel2 = @(d) exp(-d.^2/dscale(i));
    rho_bar = sum(corr_kernel2(dMids).*nrnCount)/sum(nrnCount');
    tt(i,:) = SIG0*SIG_N(t*rho_bar);
end

idcs = interp1(dscale,1:length(dscale),[5,13],'nearest');

figureNB;
    clrs = clrsPT.sequential(6);
    clrs = clrs(3:end,:);
    for i = 1:length(idcs)
        plot(t,tt(idcs(i),:),'LineWidth',1,'color',clrs(i,:));
        hold on;
    end
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    line(get(gca,'xlim'),[A1,A1],'color','k','LineStyle','--','LineWidth',1);
    line(get(gca,'xlim'),[A2,A2],'color','k','LineStyle','--','LineWidth',1);
    ylabel(['Total EEG power (' char(956) 'V^2)'])
    xlabel('\rho_{max}')
    rho1 = interp1(tt(idcs(1),:),t,A2)
    rho2 = interp1(tt(idcs(2),:),t,A2)
    plot(rho1,A2,'.k','MarkerSize',20);
    plot(rho2,A2,'.k','MarkerSize',20);