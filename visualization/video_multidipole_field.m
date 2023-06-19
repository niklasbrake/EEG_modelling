load('E:\Research_Projects\004_Propofol\presentations\_resources\multidiple_field\electric_field.mat')
load('E:\Research_Projects\004_Propofol\manuscript\Version3\Data\cortical_column_Hagen\morphology_segmentations.mat')
mData = nrnSegs.L23E_oi24rpy1;
figureNB(4.1,4.75);
axes('Position',[0,0,1,1]);
    line(mData.x',mData.y',mData.z','color','k');
    xlim([-200,200]);
    ylim([-200,400]);
    gcaformat
    set(gca,'DataAspectRatio',[1,1,1]);
    hold on;
    S = mesh(X,0*Y,Y,0*Y,'FaceAlpha',0.5,'FaceColor','interp','LineStyle','none');
    view([0,0])
    set(gca,'CLim',[-1e-2,1e-2]);
    % colorbar('color','w');
    % colormap(clrsPT.diverging(100));
    colormap(clrsPT.diverging(100))
    i= 5e3;
    S.CData = reshape(d(i,:),[M,M]);

v = VideoWriter('E:\Research_Projects\004_Propofol\presentations\_resources\multidiple_field\multidipole_field.avi','Motion JPEG AVI');
v.FrameRate = 60;
v.Quality = 30;
open(v);
frame = getframe(gcf);
for i = 5001:2:8600
    S.CData = reshape(d(i,:),[M,M]);
    drawnow;
    frame = getframe(gcf);
    writeVideo(v,frame);
    % title(sprintf('Time = %.1f ms',time(i)),'color','w');
end
close(v);
