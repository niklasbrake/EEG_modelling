load(fullfile(dataFolder,'simulations_synapse_dipole_orientation.mat'));

dMag = squeeze(max(vecnorm(dipoles,2,2)));
dMag./max(dMag)*3;

dipole_moment = squeeze(nanmedian(dipoles./vecnorm(dipoles,2,2)))';
dipole_moment = [dipole_moment(:,1),dipole_moment(:,3),-dipole_moment(:,2)];
dipole_moment(find(E_or_I==0),:) = -dipole_moment(find(E_or_I==0),:);

clrs = [230, 25, 75;
    60, 180, 75];

[theta0,phi0] = cart2sph(synapse_position(:,1),synapse_position(:,2),synapse_position(:,3));
[theta1,phi1] = cart2sph(dipole_moment(:,1),dipole_moment(:,2),dipole_moment(:,3));
synapse_position_pol = [theta0,phi0];
PD_pol = [theta1,phi1];

% gIdcs = cellfun(@(x)isempty(strfind(x,'E')),{network.neurons.mType});
gIdcs = 1-E_or_I;
clrs = [0,0,1;1,0,0];
clrs = clrs(gIdcs+1,:);

% Remove points close to the corners of the plot (to avoid variation crossing periodic boundary conditions being displayed)
idcs1 = find(~and(PD_pol(:,1)<-0.85*pi,synapse_position_pol(:,1)>0.85*pi));
idcs2 = find(~and(PD_pol(:,2)<-0.85*pi/2,synapse_position_pol(:,2)>0.85*pi/2));
idcs = intersect(idcs1,idcs2);
% Room for label
idcs1 = find(~and(PD_pol(:,1)>0.85*pi,synapse_position_pol(:,1)<0));
idcs2 = find(~and(PD_pol(:,2)>0.75*pi/2,synapse_position_pol(:,2)<0.5));
idcs = intersect(idcs,intersect(idcs1,idcs2));
% Make data sparse for display purposes
idcs = idcs(randperm(size(idcs,1),400));

fig = figureNB(6.3,3.5);
ax(1) = axes('Position',[0.16,0.3,0.3,0.55]);
ax(2) = axes('Position',[0.6,0.3,0.3,0.55]);
for i = 1:2
    axes(ax(i));
    scatter(synapse_position_pol(idcs,i),PD_pol(idcs,i),dMag(idcs),'k','filled'); hold on;
    FT = fitlm(synapse_position_pol(idcs,i),PD_pol(idcs,i),'Weights',dMag(idcs))
end
axes(ax(1));
    gcaformat;
    title('Azimuth','FontSize',7,'FontWeight','normal');
    ylabel('Dipole moment')
    xlim([-pi,pi]); xticks([-pi,0,pi]);
    ylim([-pi,pi]); yticks([-pi,0,pi]);
    xticklabels({['-' '\pi'],'0',['\pi']});
    yticklabels({['-' '\pi'],'0',['\pi']});
    line([-pi,pi],[-pi,pi],'color','r');
axes(ax(2));
    gcaformat;
    title('Elevation','FontSize',7,'FontWeight','normal');
    xl = xlabel('Synapse location (rel. soma)');
    xl.Position = [-2.2,-2.5,-1];
    xlim([-pi,pi]/2); xticks([-pi/2,0,pi/2]);
    ylim([-pi,pi]/2); yticks([-pi/2,0,pi/2]);
    xticklabels({['-' '\pi' '/2'],'0',['\pi' '/2']})
    yticklabels({['-' '\pi' '/2'],'0',['\pi' '/2']})
    line([-pi,pi]/2,[-pi,pi]/2,'color','r');

%{

figureNB(3.5,2);
    plot(dipoles(:,:,2),'LineWidth',1,'color',[192/255,0,0]);
    axis off;
    set(gcf,'color','none');
figureNB(4.2,3.8);
    [xs,ys,zs] = sphere(40);
    mesh(xs,ys,zs,0*zs,'EdgeColor','none','FaceColor','flat','FaceLighting','gouraud','FaceAlpha',0.2)
    colormap(gray)
    hold on;
    light
    lighting gouraud
    material default

    plot3(synapse_position(2,1),synapse_position(2,2),synapse_position(2,3),'.b','MarkerSize',15)

    Q = quiver3(0,0,0,dipole_moment(2,1),dipole_moment(2,2),dipole_moment(2,3),'-r','MarkerSize',20,'LineWidth',2,'MaxHeadSize',2);
    set(gca,'DataAspectRatio',[1,1,1])
    axis off


%}
