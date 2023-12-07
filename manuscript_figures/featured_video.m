n1 = load('E:\Research_Projects\004_Propofol\data\simulations\raw\critical_dipole_correlation\example_simulation\simulation\L23E_oi24rpy1_multi_dipoles.mat')
n1.x(:,1) = n1.x(:,1)-250;

n2 = load('E:\Research_Projects\004_Propofol\data\simulations\raw\critical_dipole_correlation\example_simulation\simulation\L23I_oi38lbc1_multi_dipoles.mat')
n2.x(:,1) = n2.x(:,1)+250;

M1 = 10;
M2 = 5;
[X,Y] = meshgrid(linspace(-533.33,533.33,M1),linspace(-200,400,M2));
pts = [X(:),Y(:),0*Y(:)];

d = zeros(size(n1.dipoles,3),size(pts,1));
for k = 1:size(pts,1);
    waitbar(k/size(pts,1));
    pt = pts(k,:);

    v1 = n1.x-pt;
    r1 = vecnorm(v1,2,2);
    vn1 = v1./r1;

    v2 = n2.x-pt;
    r2 = vecnorm(v2,2,2);
    vn2 = v2./r2;

    d(:,k) = sum(squeeze(sum(n1.dipoles.*vn1,2))./r1.^2)+sum(squeeze(sum(n2.dipoles.*vn2,2))./r2.^2);
end

% load('E:\Research_Projects\004_Propofol\presentations\_resources\multidiple_field\electric_field.mat')
load('E:\Research_Projects\004_Propofol\manuscript\Version3\Data\cortical_column_Hagen\morphology_segmentations.mat')
mData1 = nrnSegs.L23E_oi24rpy1;
mData2 = nrnSegs.L23I_oi38lbc1;

v = VideoWriter('E:\Research_Projects\004_Propofol\presentations\_resources\multidiple_field\dyad_multidipole_field.avi','Motion JPEG AVI');
v.FrameRate = 60;
v.Quality = 30;
open(v);



figureNB(12,6.75);
axes('Position',[0,0,1,1]);
    line(mData1.x'-250,mData1.y',mData1.z','color','w');
    line(mData2.x'+250,mData2.y',mData2.z','color','w');
    xlim([-533.33,533.33]);
    ylim([-200,400]);
    gcaformat_dark
    view([0,0])
    set(gca,'DataAspectRatio',[1,1,1]);
    hold on;
    S = mesh(X,0*Y,Y,0*Y,'FaceAlpha',0.5,'FaceColor','interp','LineStyle','none');
    set(gca,'CLim',[-1e-3,1e-3]);
    % colorbar('color','w');
    % colormap(clrsPT.diverging(100));
    colormap(flipud(1-clrsPT.diverging(100)))

% frame = getframe(gcf);
txt = text(0,10,400,sprintf('%d',0),'Color','w','HorizontalAlignment','center','VerticalAlignment','top')
for i = 13700:5:32e3
    S.CData = reshape(d(i,:),[M2,M1]);
    txt.String = sprintf('%d',i);
    drawnow;
    % frame = getframe(gcf);
    % writeVideo(v,frame);
    % title(sprintf('Time = %.1f ms',time(i)),'color','w');
end
close(v);
