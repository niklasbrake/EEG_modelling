function figureS2(dataFolder)

    if(nargin<1)
        error('Path to data required as input argument. Data can be downloaded from link in README file.');
    end


    load(fullfile(dataFolder,'simulations','synthetic_spikes','simulations_synapse_dipole_orientation.mat'));

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

    gIdcs = 1-E_or_I;
    clrs = [0,0,1;1,0,0];
    clrs = clrs(gIdcs+1,:);

    mIDs = findgroups(mTypes);
    mTypes = unique(mTypes);

    for j = 1:11
        fig = figureNB(6,3);
        for i = 1:2
            subplot(1,2,i);
            idcs = find(mIDs==j);
            scatter(synapse_position_pol(idcs,i),PD_pol(idcs,i),dMag(idcs),'k','filled'); hold on;
            FT = fitlm(synapse_position_pol(idcs,i),PD_pol(idcs,i),'Weights',dMag(idcs))
        end
        subplot(1,2,1);
            gcaformat;
            xlim([-pi,pi]); xticks([-pi,0,pi]);
            ylim([-pi,pi]); yticks([-pi,0,pi]);
            xticklabels({sprintf('%c%c',char(45),char(960)),'0',sprintf('%c',char(960))})
            yticklabels({sprintf('%c%c',char(45),char(960)),'0',sprintf('%c',char(960))})
            line([-pi,pi],[-pi,pi],'color','r');
            text(pi,pi+0.3,strrep(mTypes{j},'_','\_'),'VerticalAlignment','bottom','FontSize',7,'FontWeight','bold')
            axis square;
        subplot(1,2,2);
            gcaformat;
            xlim([-pi,pi]/2); xticks([-pi/2,0,pi/2]);
            ylim([-pi,pi]/2); yticks([-pi/2,0,pi/2]);
            xticklabels({sprintf('%c%c/2',char(45),char(960)),'0',sprintf('%c/2',char(960))})
            yticklabels({sprintf('%c%c/2',char(45),char(960)),'0',sprintf('%c/2',char(960))})
            line([-pi,pi]/2,[-pi,pi]/2,'color','r');
            axis square;
    end
end