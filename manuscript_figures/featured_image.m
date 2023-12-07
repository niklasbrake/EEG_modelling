function featured_image(n1,n2,i0)

load('E:\Research_Projects\004_Propofol\manuscript\Version3\Data\cortical_column_Hagen\morphology_segmentations.mat')

M = 100;
[X,Y] = meshgrid(linspace(-533.33,533.33,M),linspace(-200,400,M));
pts = [X(:),Y(:),0*Y(:)];

% i0 = 11190
if(nargin<3)
    i0 = 13990
end


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

mData1 = nrnSegs.L23E_oi24rpy1;
mData2 = nrnSegs.L23I_oi38lbc1;
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
    S.CData = reshape(d,[M,M]);
    set(gca,'CLim',[-1e-2,1e-2]);
    % colorbar('color','w');
    % colormap(clrsPT.diverging(100));
    colormap(flipud(1-clrsPT.diverging(100)))
