n1 = load('E:\Research_Projects\004_Propofol\data\simulations\raw\example_embedding\simulation\L23E_oi24rpy1_multi_dipoles.mat')
n1.x(:,1) = n1.x(:,1)-250;

n2 = load('E:\Research_Projects\004_Propofol\data\simulations\raw\example_embedding\simulation\L23I_oi38lbc1_multi_dipoles.mat')
n2.x(:,1) = n2.x(:,1)+250;


load('E:\Research_Projects\004_Propofol\manuscript\Version3\Data\cortical_column_Hagen\morphology_segmentations.mat')

M1 = 1200;
M2 = 675;
[X,Y] = meshgrid(linspace(-533.33,533.33,M1),linspace(-200,400,M2));
pts = [X(:),Y(:),0*Y(:)];

i0 = 25256;

d = zeros(1,size(pts,1));
for k = 1:size(pts,1);
    waitbar(k/size(pts,1));
    pt = pts(k,:);
    
    v1 = n1.x-pt;
    r1 = vecnorm(v1,2,2);
    vn1 = v1./r1;
    
    v2 = n2.x-pt;
    r2 = vecnorm(v2,2,2);
    vn2 = v2./r2;

    d(k) = sum(squeeze(sum(n1.dipoles(:,:,i0).*vn1,2))./r1.^2)+sum(squeeze(sum(n2.dipoles(:,:,i0).*vn2,2))./r2.^2);
end
CM = [1,1,1;flipud(1-clrsPT.diverging(100))];
mData1 = nrnSegs.L23E_oi24rpy1;
mData2 = nrnSegs.L23I_oi38lbc1;
d = min(max(d,-1e-2),1e-2);

figureNB(12,6.75);
axes('Position',[0,0,1,1]);
    gcaformat_dark
    colormap(CM);
    set(gca,'CLim',[-1.1e-2,1.1e-2]);
    hold on;
    render_neuron_morphology(1,-250);
    render_neuron_morphology(2,250);
    xlim([-533.33,533.33]);
    ylim([-200,400]);
    gcaformat_dark
    view([0,0])
    set(gca,'DataAspectRatio',[1,1,1]);
    hold on;
    S = mesh(X,0*Y,Y,0*Y,'FaceAlpha',0.5,'FaceColor','interp','LineStyle','none');
    S.CData = reshape(d,[M2,M1]);

    print(gcf,'E:\Research_Projects\004_Propofol\manuscript\Version3\Figures\png\featured_image.bmp','-dbmp','-r600');

function fig = render_neuron_morphology(mType,offset)

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
    data(:,3:5) = data(:,3:5)-mean(xSoma(:,1:3));
    xSoma(:,1:3) = xSoma(:,1:3)-mean(xSoma(:,1:3));

    data(:,3) = data(:,3)+offset;
    xSoma(:,1) = xSoma(:,1)+offset;

    gridSize = @(d) floor(d.^2./(3+d.^2)*18+3);
    % fig = figureNB(12,18);
    % axes('Position',[0,0,1,1]);
    view([0,0]);
    hold on
    set(gca,'DataAspectRatio',[1,1,1]);
    gcaformat_dark
    axis off;

    xSynapses = cell(length(segs)+1,1);
    ei = cell(length(segs)+1,1);
    for i = 1:length(segs)
        idcs = segs{i};
        if(max(data(idcs,6))>1)
            xSegment = data(idcs,3:6);
        else
            xSegment = data(idcs(1),3:6);
            pR = data(idcs(1),6);
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
        N = gridSize(max(xSegment(:,end)));
        [xSynapses{i},ei{i}] = render_segment(xSegment,N);
        drawnow;
    end
    N = gridSize(max(xSoma(:,end)));
    [xSynapses{end+1},ei{end+1}] = render_segment(xSoma,N);
    xSynapses = cat(1,xSynapses{:});
    ei = cat(1,ei{:});
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

    K = -1.1e-2+0*exp(-Yall.^2/80^2);
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
end