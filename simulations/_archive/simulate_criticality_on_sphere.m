tmax = 2;

branchNo = 0.98;
N = 2e3;
K1 = linspace(0,1,N);
elevation = asin(2*rand(N,1)-1);
azimuth = 2*pi*rand(N,1);
[x,y,z] = sph2cart(azimuth,elevation,1);

dt = 4e-3;
nNeigh = 8;
pmax = branchNo/nNeigh;
C = zeros(N*nNeigh,4);
for i = 1:N
    network_simulation.haversine_distance2([elevation(i),azimuth(i)],[elevation,azimuth]);
    [~,I] = sort(network_simulation.haversine_distance2([elevation(i),azimuth(i)],[elevation,azimuth]),'ascend');
    C(nNeigh*(i-1)+1:nNeigh*i,1) = i*ones(nNeigh,1);
    C(nNeigh*(i-1)+1:nNeigh*i,2) = I(2:nNeigh+1);
    C(nNeigh*(i-1)+1:nNeigh*i,3) = rand(nNeigh,1)*pmax*2;
    C(nNeigh*(i-1)+1:nNeigh*i,4) = dt*2*rand(nNeigh,1); % Random transmission speed
end

[ids,ts,ei,C2] = simulatespikes_det(N,branchNo,tmax,C,nNeigh);


figureNB;
for i = 1:4
    subplot(2,2,i);
    histogram(C(:,i));
    hold on;
    histogram(C2(:,i));
end


filename = 'testnew51.gif';
fig = figureNB;
S = scatter3(x,y,z,ones(N,1),zeros(N,1),'filled');
set(gca,'CLim',[0,1]);
set(gca,'DataAspectRatio',[1,1,1]);
colormap([0.3,0.3,0.3;1,1,1]);
gcaformat_dark;
axis off;
axis vis3d
for t = 0:1:tmax*1e3/2
    idcs = ids(find(and(ts>t,ts<t+dt*1e3)));
    S.CData(idcs) = 1;
    S.SizeData(idcs) = 5;
    camorbit(1,0,'coordsys',[0,0,1])
    drawnow;
    frame = getframe(fig);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if t == 0;
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',1/30);
    end
    S.CData(idcs) = 0;
    S.SizeData(idcs) = 1;
end
