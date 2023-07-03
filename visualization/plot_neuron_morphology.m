load('E:\Research_Projects\004_Propofol\manuscript\Version3\Data\cortical_column_Hagen\matlab_morph_data\L23E_connections.mat')
connections = connections+1;
load('E:\Research_Projects\004_Propofol\manuscript\Version3\Data\cortical_column_Hagen\matlab_morph_data\L23E.mat')
X = csvread('E:\Research_Projects\004_Propofol\manuscript\Version3\Data\cortical_column_Hagen\matlab_morph_data\L23E.csv');
for i=  1:size(X,1)
    segs{i} = X(i,X(i,:)~=0);
end


data(:,6) = max(data(:,6),0.5);
xSoma = data(data(:,2)==1,3:6);

gridSize = @(d) floor(d.^2./(3+d.^2)*18+3);
figureNB(12,18);
axes('Position',[0,0,1,1]);
view([0,0]);
hold on
set(gca,'DataAspectRatio',[1,1,1]);
axis off;
set(gca,'CLim',[0,1]);
colormap('gray');
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

    % xSegment = [x0;xSegment];

    N = gridSize(max(xSegment(:,end)));
    render_segment(xSegment,N)
    drawnow;
end
N = gridSize(max(xSoma(:,end)));
render_segment(xSoma,N);


function render_segment(xSegment,N)
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

    K = 1-exp(-Yall.^2/80^2);
    surf(Xall,Yall,Zall,K,'LineStyle','none','FaceLighting','gouraud');
    % surf(Xall,Yall,Zall,'LineStyle','none');
end