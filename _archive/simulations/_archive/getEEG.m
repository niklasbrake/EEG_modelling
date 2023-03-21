function [W,idx] = getEEG(Q,sa,idx)
    if(nargin<2)
        load('E:\Research_Projects\004_Propofol\Modelling\head_models\sa_nyhead.mat')
    end
    if(nargin<3)
        idx = randi(size(sa.cortex75K.vc,1));
    end
    % Orient dipole "up" direction to be normal to cortical surface
    vz = sa.cortex75K.normals(idx,:);
    [vx,vy] = getOrthBasis(vz);
    % Lead field for Cz recording
    L0 = squeeze(sa.cortex75K.V_fem(49,idx,:))'; % czIDX = 49
    q = (Q(:,1).*vx+Q(:,2).*vy+Q(:,3).*vz)*1e-12; % nAum -> mAm
    W = sum(L0.*q,2)*1e6; % V -> uV
end



