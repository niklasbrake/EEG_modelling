function figureS1(dataFolder)

    if(nargin<1)
        error('Path to data required as input argument. Data can be downloaded from link in README file.');
    end

    for i = 1:11
        render_neuron_morphology(i,dataFolder);
    end

end

function fig = render_neuron_morphology(mType,dataFolder)

    mTypes ={'L23E_oi24rpy1';'L23I_oi38lbc1';'L4E_53rpy1';'L4E_j7_L4stellate';'L4I_oi26rbc1';'L5E_j4a';'L5E_oi15rpy4';'L5I_oi15rbc1';'L6E_51_2a_CNG';'L6E_oi15rpy4';'L6I_oi15rbc1'};

    if(isnumeric(mType))
        mType = mTypes{mType}
    end

    folder = fullfile(dataFolder,'cortical_column_Hagen','matlab_morph_data');

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
end