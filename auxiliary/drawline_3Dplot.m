function L = drawline_3Dplot(x0,thet,R,H)
    VW = get(gca,'view');
    [a1,a2,a3] = sph2cart((VW(1)-90)*pi/180,VW(2)/180*pi,1);
    e0 = [a1,a2,a3];
    e2 = camup;
    e1 = cross(e0,-e2);
    L = [zeros(1,3);R*(e1*cos(thet)+e2*sin(thet))]'+x0(:);
    if(nargin==4)
        L = L + e0'*H;
    end
    L = line(L(1,:),L(2,:),L(3,:));
