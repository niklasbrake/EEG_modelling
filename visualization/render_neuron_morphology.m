function fig = render_neuron_morphology(mType)

mTypes ={'L23E_oi24rpy1';'L23I_oi38lbc1';'L4E_53rpy1';'L4E_j7_L4stellate';'L4I_oi26rbc1';'L5E_j4a';'L5E_oi15rpy4';'L5I_oi15rbc1';'L6E_51_2a_CNG';'L6E_oi15rpy4';'L6I_oi15rbc1'};

if(isnumeric(mType))
    mType = mTypes{mType}
end

folder = 'E:\Research_Projects\004_Propofol\manuscript\Version3\Data\cortical_column_Hagen\matlab_morph_data';

load(fullfile(folder,[mType '_connections.mat']));
connections = connections+1;

X = csvread(fullfile(folder,[mType '_segments.csv']));
for i=  1:size(X,1)
    segs{i} = X(i,X(i,:)~=0);
end



data(:,6) = max(data(:,6),0.5);
xSoma = data(data(:,2)==1,3:6);
data(:,3:5) = data(:,3:5) - mean(xSoma(:,1:3));
xSoma(:,1:3) = xSoma(:,1:3)-mean(xSoma(:,1:3));

gridSize = @(d) floor(d.^2./(3+d.^2)*18+3);
fig = figureNB(12,18);
axes('Position',[0,0,1,1]);
view([0,0]);
view([0,0]);
hold on
set(gca,'DataAspectRatio',[1,1,1]);
% gcaformat_dark
axis off;
set(gca,'CLim',[0,1]);
colormap(flipud(gray(1e3)));

xSynapses = cell(length(segs)+1,1);
ei = cell(length(segs)+1,1);
for i = 1:length(segs)
    idcs = segs{i};
    % M(i) = max(data(idcs,6));
    if(max(data(idcs,6))>1)
        xSegment = data(idcs,3:6);
    else
        xSegment = data(idcs(1),3:6);
        pR = data(idcs(1),6);
        % if(max(data(idcs,6))1)
        %     for j = 2:length(idcs)-1
        %         if(abs(data(idcs(j),6)-pR)./pR>0.15)
        %             pR = data(idcs(j),6)-pR;
        %             xSegment = [xSegment;data(idcs(j),3:6)];
        %         end
        %     end
        % end
        xSegment = [xSegment;data(idcs(end),3:6)];
    end


    j0 = find(connections(:,2)==i);
    i0 = connections(j0,1);
    if(i0==0)
        x0 = xSoma(1,:);
    else
        x0 = data(segs{i0},3:6);
        x0 = x0(end-1,:);
    end

    % xSegment = [x0;xSegment];

    N = gridSize(max(xSegment(:,end)));
    [xSynapses{i},ei{i}] = render_segment(xSegment,N);
    drawnow;
end
N = gridSize(max(xSoma(:,end)));
[xSynapses{end+1},ei{end+1}] = render_segment(xSoma,N);
xSynapses = cat(1,xSynapses{:});
ei = cat(1,ei{:});


zl = get(gca,'zlim');
xl = get(gca,'xlim');
if(xl(1)>-100)
    xl(1) = -100;
end
if(xl(2)<100)
    xl(2) = 100;
end
zl(1) = zl(1)-20;
xlim(xl);
zlim(zl);
line([-100,100],[0,0],[1,1]*zl(1),'LineWidth',2,'color','k');
print(fig,fullfile('E:\Research_Projects\004_Propofol\presentations\_resources\neuron_morphologies',[mType '.svg']),'-dsvg')
% close(fig);
return;
set(gca,'CLim',[-1,1])
colormap([1,0,0;0,0,0;0,0,1]);
idcs = randsample(size(xSynapses,1),1e3);
S = scatter3(xSynapses(idcs,1),xSynapses(idcs,2),xSynapses(idcs,3),0.5,ei(idcs)*2-1,'filled');


% To get the section -> segment map, go to python and do
% for i in range(cell.totnsegs):
%    secNames.append(cell.get_idx_name(i)[1])
% Then average over each segment to get the section voltages
load('E:\Research_Projects\004_Propofol\data\simulations\raw\dendrite_voltages\LFPy\voltage_data.mat')
V = [vDend,vSoma];
V2 = resample(V,1e3,16e3);
t = (1:size(V2,1));

ch = flipud(get(gca,'Children'));
for i = 1:length(ch)
    alpha{i} = 1+0*ch(i).CData;
end
% F = (V2 - -65)/(-55--65);
% F = min(max(F,0),1);

set(gca,'CLim',[-70,-50]);

CM = jet(1001);
colormap(flipud([CM(502:end,:);CM(1:500,:)]).*(1-exp(-linspace(-5,5,1e3).^2)'));

v = VideoWriter('newfile.avi','Motion JPEG AVI');
v.Quality = 95;
open(v)
for i = 800:1000
    for j = 1:length(ch)
        ch(j).CData = max(alpha{j}*V2(i,j),-62);
    end
    title(sprintf('%.1f ms',t(i)),'color','w');
    drawnow;
   frame = getframe(gcf);
   writeVideo(v,frame);
end
close(v);

end
function [xSynapses,ei] = render_segment(xSegment,N)
    Xall = [];
    Yall = [];
    Zall = [];
    for i = 1:size(xSegment,1)-1
        x0 = xSegment(i,1:3);
        [X,Y,Z] = cylinder(xSegment(i:i+1,4),N);
        vz = diff(xSegment(i:i+1,1:3));
        a = norm(diff(xSegment(i:i+1,1:3)));
        [vx,vy] = getOrthBasis(vz/a);

        X2 = X;
        Y2 = Y;
        Z2 = Z;
        for j = 1:size(X,2)
            temp = X(:,j).*vx + Y(:,j).*vy + Z(:,j).*vz;
            X2(:,j) = temp(:,1);
            Y2(:,j) = temp(:,2);
            Z2(:,j) = temp(:,3);
        end

        Xall = [Xall;X2+x0(1)];
        Yall = [Yall;Y2+x0(2)];
        Zall = [Zall;Z2+x0(3)];
    end
    x0 = xSegment(end,1:3);
    % Xall = [Xall;X2+x0(1)];
    % Yall = [Yall;Y2+x0(2)];
    % Zall = [Zall;Z2+x0(3)];

    K = 1+0*exp(-Yall.^2/80^2);
    surf(Xall,Yall,Zall,K,'LineStyle','none','FaceLighting','gouraud');

    xSynapses= [];
    ei= [];
    return
    d = cumsum([0;vecnorm(diff([Xall(:,1),Yall(:,1),Zall(:,1)]),2,2)]);
    % Xall

    nE = poissrnd(d(end));
    idcs = rand(nE,1);
    xE = interp1(d/max(d),Xall,idcs);
    yE = interp1(d/max(d),Yall,idcs);
    zE = interp1(d/max(d),Zall,idcs);
    idcs = 2*pi*rand(nE,1);
    xE2 = nan(length(idcs),1);
    yE2 = nan(length(idcs),1);
    zE2 = nan(length(idcs),1);
    for j = 1:length(idcs)
        xE2(j) = interp1(linspace(0,2*pi,N+1),xE(j,:),idcs(j));
        yE2(j) = interp1(linspace(0,2*pi,N+1),yE(j,:),idcs(j));
        zE2(j) = interp1(linspace(0,2*pi,N+1),zE(j,:),idcs(j));
    end


    nI = poissrnd(0.15*d(end));
    idcs = rand(nI,1);
    xI = interp1(d/max(d),Xall,idcs);
    yI = interp1(d/max(d),Yall,idcs);
    zI = interp1(d/max(d),Zall,idcs);
    idcs = 2*pi*rand(nI,1);
    xI2 = nan(length(idcs),1);
    yI2 = nan(length(idcs),1);
    zI2 = nan(length(idcs),1);
    for j = 1:length(idcs)
        xI2(j) = interp1(linspace(0,2*pi,N+1),xI(j,:),idcs(j));
        yI2(j) = interp1(linspace(0,2*pi,N+1),yI(j,:),idcs(j));
        zI2(j) = interp1(linspace(0,2*pi,N+1),zI(j,:),idcs(j));
    end

    xSynapses = [xE2(:),yE2(:),zE2(:); ...
                xI2(:),yI2(:),zI2(:)];
    ei = [zeros(nE,1);ones(nI,1)];

    % plot3(Xall(iI),Yall(iI),Zall(iI),'.b','MarkerSize',5);

    % surf(Xall,Yall,Zall,'LineStyle','none');
end