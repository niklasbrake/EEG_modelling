% load('E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\raw\time_series_all_channels.mat')
y = detrend(TimeDomainAligned(40e4+1:41e4,2,1));
N = length(y);
fr = 1024/N;
fs = [0:fr:1024/2 [1024/2+fr:fr:1024-fr]-1024];

yk = y(:)'*exp(-sqrt(-1)*2*pi*((1:N)'-1).*((1:N)-1)/N);

yi = zeros(N,N);
for i = 1:N
	yi(:,i) = 1/N*yk(i).*exp(sqrt(-1)*2*pi*(i-1)*((1:N)'-1)/N);
end

for i = 1:N
	amp(i) = abs(yk(i))/N;
	phase(i) = acos(real(yk(i))/abs(yk(i)));
	sgn(i) = sign(imag(yk(i)));
end

% Get most important frequencies
F2 = 55;

[a,b] = sort(amp(1:floor(F2/fr)),'descend');
b = b(1:12);
b = sort(b);

y_approx = real(yi(:,b));
y2 = y/range(y)*range(sum(y_approx,2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = 'E:\Research_Projects\004_Propofol\presentations\2023_HBHL\expalinPSD.gif';

clrs(1,:) = [255,212,0]/255;
clrs(2,:) = [0,176,240]/255;
clrs(3,:) = [254,158,68]/255;
clrs(4,:) = [1,194,154]/255;
clrs(5,:) = [97,113,217]/255;
clrs(6,:) = [54,199,208]/255;
clrs(7,:) = [253,99,65]/255;


v = VideoWriter('newfile.avi','Motion JPEG AVI');
v.Quality = 95;
open(v);

fig = figure('color','k','position',[705.8000  317.0000  656.4000  381.2000]);
h = axes('Position',[0.05,0.02,0.91,0.96]); hold on;
set(gca,'colororder',clrs);
iValues = linspace(0,1,100);
iValues2 = linspace(0,1,25);
ylim2 = N/1024;

plot3(0*y2,[1:N]/1024,y2); hold on;
F = fill3([0,0,0,0],[0,ylim2,ylim2,0],[-10,-10,10,10],'w','EdgeColor','w');
F.FaceAlpha = 0.15;
xlim([0,1]);
zlim([-15,15]);
ylim([0,ylim2]);
view([-90,0])
set(gca,'color','k');
set(get(gca,'xaxis'),'visible','off')
set(get(gca,'yaxis'),'visible','off')
set(get(gca,'zaxis'),'visible','off')

% Capture the plot as an image 
drawnow;
im = frame2im(getframe(fig)); 
[imind,cm] = rgb2ind(im,256); 
imwrite(imind,cm,filename,'gif', 'Loopcount',0,'DelayTime',3);

pl0 = plot3(0*y,[1:N]/1024,sum(y_approx,2));
drawnow;
im = frame2im(getframe(fig)); 
[imind,cm] = rgb2ind(im,256); 
imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.5);

v0 = [-90,0];
v1 = [-37.5,30];
for i = iValues2
	view(v0+i*(v1-v0));
	drawnow;
	im = frame2im(getframe(fig)); 
	[imind,cm] = rgb2ind(im,256); 
	imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',1/25);
end

delete(pl0);
for i = 1:length(b)
	pl(i) = plot3(0*y,[1:N]/1024,y_approx(:,i));
end
for i = iValues2
	for j = 1:length(b)
		pl(j).XData = i*(j/length(b)) + 0*y;
	end
	drawnow;
	im = frame2im(getframe(fig)); 
	[imind,cm] = rgb2ind(im,256); 
	imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',1/25);
end
imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',5);

F = fill3([0,1,1,0],[0,0,0,0],[-10,-10,10,10],'w','EdgeColor','w');
F.FaceAlpha = 0;
Lxax = line([0,1],[0,0],[0,0],'Color','w');
for i = 1:length(b)
	L(i) = line([i/length(b) i/length(b)],[0,0],[0,amp(b(i))],'Color','w');
end

v0 = [-37.5,30];
v1 = [0,0];
for i = iValues2
	view(v0+i*(v1-v0));
	F.FaceAlpha = i;
	F.FaceColor = [0,0,0] + 0.15*i;
	drawnow;
	im = frame2im(getframe(fig)); 
	[imind,cm] = rgb2ind(im,256); 
	imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',1/25);
end
imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.5);


z0 = [];
for i = iValues2
	for j = 1:length(L)
		L(j).XData = [0,0] + j/length(b) + i*((b(j)*fr)/F2 - j/length(b));
		M(j) = max(L(j).ZData);
		z0(:,j) = L(j).ZData;
	end
	drawnow;
	im = frame2im(getframe(fig)); 
	[imind,cm] = rgb2ind(im,256); 
	imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',1/25);
end

iVec = [1 10:10:F2];
for j = 1:length(iVec)
	i = iVec(j);
	LM(j) = line([i/F2,i/F2],[0,0],[-0.5,0],'Color','w');
	LT(j) = text(i/F2,0,-0.5,int2str(i),'HorizontalAlignment','center', ...
		'VerticalAlignment','top','FontSize',12,'color','w');
end
LXL = text([0.5,0.5],[0,0],[-1.75,-1.75],'Frequency (Hz)','HorizontalAlignment','center', ...
		'VerticalAlignment','top','FontSize',12,'color','w');
T = text([-0.05,-0.05],[0,0],[5,5],'Power','HorizontalAlignment','center', ...
		'VerticalAlignment','top','Rotation',90,'FontSize',12,'color','w');
drawnow;
im = frame2im(getframe(fig)); 
[imind,cm] = rgb2ind(im,256); 
imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.5);

M = max(M);
for i = iValues2
	for j = 1:length(L)
		L(j).ZData = z0(:,j) * (1+i*(9.5/M-1));
	end
	drawnow;
	im = frame2im(getframe(fig)); 
	[imind,cm] = rgb2ind(im,256); 
	imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',1/25);
end

xValues = [fr:fr:F2]/F2;
z0 = amp(2:length(xValues)+1).^2;
z0 = 9.5*(z0-min(z0))/(max(z0)-min(z0));
p=plot3(xValues,0*amp(2:length(xValues)+1),z0,'w');
delete(L);
posXLab = get(LXL,'Position');
posXLab = posXLab{1};
posYLab = get(T,'Position');
posYLab = posYLab{1};
y0 = p.ZData;
for i = iValues2
	y0 =  (9.5+i*10)*(y0-min(y0))/(max(y0)-min(y0))-i*10;
	p.ZData = y0;
	Lxax.ZData = [0,0] + i*(-10);
	posXLab(3) = -1.75 - i*10;
	posYLab(3) = 5 - i*5;
	set(LXL,'Position',posXLab);
	set(T,'Position',posYLab);
	for j =1:length(LM)
		LT(j).Position(3) = LT(j).Position(3) - iValues2(2)*10*(i>0);
		LM(j).ZData = [-0.5,0] - i*10;
	end
	drawnow;
	im = frame2im(getframe(fig)); 
	[imind,cm] = rgb2ind(im,256); 
	imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',1/25);
end

y0 = p.ZData;
y1 = log(smoothdata(amp(2:length(xValues)+1).^2,'movmean',3));
y0 = 19.5*(y0-min(y0))/(max(y0)-min(y0))-10;
y1 = 19.5*(y1-min(y1))/(max(y1)-min(y1))-10;
for i = iValues2
	p.ZData = y0 + i*(y1-y0);
	drawnow;
	im = frame2im(getframe(fig)); 
	[imind,cm] = rgb2ind(im,256); 
	imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',1/25);
end
set(T,'String','Log Power');

im = frame2im(getframe(fig)); 
[imind,cm] = rgb2ind(im,256); 
imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',10);
