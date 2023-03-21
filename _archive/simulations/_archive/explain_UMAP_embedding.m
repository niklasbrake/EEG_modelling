load('E:\Research_Projects\004_Propofol\Modelling\neuron_simulations\data\cortical_column_Hagen\segment_areas.mat');
mData = nrnSegs.L23E_oi24rpy1;


R = max(max(sqrt(mData.x.^2+mData.y.^2+mData.z.^2)));


pos = [mean(mData.x,2),mean(mData.y,2),mean(mData.z,2)];
posN = R*pos./vecnorm(pos,2,2);

i = randi(size(mData.x,1));
figureNB(17,17);
    [X,Y,Z] = sphere(20);
    M = mesh(R*X,R*Y,R*Z);
    M.EdgeColor = [0.3,0.3,0.3];
    M.FaceColor = [0.4,0.4,0.4];
    M.FaceAlpha = 0.2;

    hold on;
    line(mData.x',mData.y',mData.z','color','w');
    gcaformat_dark;
    axis off;
    set(gca,'DataAspectRatio',[1,1,1])
    scatter3(pos(i,1),pos(i,2),pos(i,3),15,'y','filled');
    scatter3(posN(i,1),posN(i,2),posN(i,3),15,'y','filled');
    line([0,posN(1)],[0,posN(2)],[0,posN(3)],'color','y');


figureNB(17,17);
    [X,Y,Z] = sphere(20);
    M = mesh(R*X,R*Y,R*Z);
    M.EdgeColor = [0.3,0.3,0.3];
    M.FaceColor = [0.4,0.4,0.4];
    M.FaceAlpha = 0.2;

    hold on;
    scatter3(posN(:,1),posN(:,2),posN(:,3),5,'w','filled');
    gcaformat_dark;
    axis off;
    set(gca,'DataAspectRatio',[1,1,1])




for S = [0,0.05,0.1,0.2,0.5,1]
fig = figureNB(10,10);
axes('Position',[0,0,1,1])
    [X,Y,Z] = sphere(25);
    M = mesh(X,Y,Z);
    M.LineStyle = ':';
    M.EdgeColor = [0.4,0.4,0.4];
    M.FaceColor = [0.4,0.4,0.4];
    M.FaceAlpha = 0.2;
    hold on;
    gcaformat_dark;
    axis off;
    set(gca,'DataAspectRatio',[1,1,1])

    X = repmat([0.1282,-0.671,0.7302],[1e2,1]);
    [theta0,phi0] = cart2sph(X(:,1),X(:,2),X(:,3));
    phi0 = phi0+pi/2;
    sz = size(phi0);
    for k = 1:40
        B = pi*rand*ones(sz);
        u = 1-2*S*rand;
        idcs = (rand<0.5);
        if(idcs)
            d = acos(u);
            d = linspace(0,d,size(X,1))';
        else
            d = 2*pi-acos(u);
            d = linspace(2*pi,d,size(X,1))';
        end
        phi2 = acos(cos(d).*cos(phi0)+sin(d).*sin(phi0).*cos(B));
        A = acos((cos(d)-cos(phi2).*cos(phi0))./(sin(phi2).*sin(phi0)));
        A(d>pi) = 2*pi-A(d>pi);
        theta2 = theta0+real(A);
        [x,y,z] = sph2cart(theta2,phi2-pi/2,1);
        Xnew = [x,y,z];

        scatter3(Xnew(1,1),Xnew(1,2),Xnew(1,3),15,'w','filled');
        plot3(Xnew(:,1),Xnew(:,2),Xnew(:,3),'color',[0.6,0.6,0.6],'LineWidth',1);
        XX(k,:) = Xnew(end,:);
    end
    scatter3(XX(:,1),XX(:,2),XX(:,3),15,'y','filled');
fig.InvertHardcopy = 'off';
saveas(fig,['C:/users/brake/Desktop/s=' num2str(S) '.png']);
close(fig);
end