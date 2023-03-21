function [A,A0] = calculate_intersection_area(tri,R)

% Entire triangle is within sphere
if(all(vecnorm(tri,2,2)<R))
    A = 1/2*norm(cross(tri(1,:)-tri(2,:),tri(3,:)-tri(1,:)));
    return;
end

% Get plane of the triangle
params = null([tri,ones(3,1)])';
N = params(1:3);
rho = -params(4)/norm(N,2);

% Entire triangle is outside sphere  (type 1)
if(abs(rho)>R)
    A = 0;
    return;
end

% Calculate circle of the sliced sphere
R0 = sqrt(R^2-rho^2);
c1 = rho*N/norm(N,2);

% Calculate intersection points
pairs = [1,2;2,3;3,1];
A = sum(tri(pairs(:,1),:).^2,2) - R^2;
C = sum((tri(pairs(:,1),:)-tri(pairs(:,2),:)).^2,2);
B = sum(tri(pairs(:,2),:).^2,2) - A - C - R^2;
in = zeros(6,3);
T = zeros(6,3);
for k = 1:3
    t = roots([C(k),B(k),A(k)]);
    in(2*k-1:2*k,:) = [tri(pairs(k,1),:).*(1-t) + tri(pairs(k,2),:).*t];
    if(sum(t<0) == 2 || sum(t>1) == 2)
        % The line intersects, but segment doesn't
        % Combine this case with imaginary roots
        t = t*sqrt(-1);
    end
    T(2*k-1:2*k,:) = [pairs(k,:),t(1);pairs(k,:),t(2)];
end

% No intersection points
if(all(imag(T(:,3))~=0))
    if(norm(mean(tri),2)<R)
        % Entire sphere is inside triangle
        A = pi*R0^2;
    else
        % Entire triangle is outside sphere (type 2)
        A = 0;
    end
    return
end

% Compute vertices that define intersected area
in(T(:,3)<0,:) = tri(T(T(:,3)<0,1),:);
in(T(:,3)>1,:) = tri(T(T(:,3)>1,2),:);
in(imag(T(:,3))~=0,:) = [];
T(imag(T(:,3))~=0,:) = [];
[T,idcs] = sortrows(T);
in = in(idcs,:);

% Change coordinates to intersecting plane
[x1,x2] = getOrthBasis(params(1:3));
in_proj = [sum((in-c1).*x1,2),sum((in-c1).*x2,2)];
in_proj(end+1,:) = in_proj(1,:);

% Get area of polygon
A = sum((in_proj(2:end,1)+in_proj(1:end-1,1)).*(diff(in_proj(:,2))))/2;

% Get sign of area to determine direction of circular segments
clockwise = (A<0);

% Convert to polar coordinates
theta = mod(atan2(in_proj(:,2),in_proj(:,1)),2*pi);

% Angles that define start and end points of each circular segments
arcTheta = [];
for k = 2:2:length(theta)
    if(theta(k)~=theta(k+1))
        if(clockwise)
            if(theta(k+1)>theta(k))
                theta(k+1) = theta(k+1)-2*pi;
            end
        else
            if(theta(k+1)<theta(k))
                theta(k+1) = theta(k+1)+2*pi;
            end
        end
        arcTheta(end+1,:) = theta(k:k+1);
    end
end
% If entire intersection is one circular segment
if(size(arcTheta,1)==1)
    arcTheta = mod(arcTheta,2*pi);
    alph1 = mod(0.5*(arcTheta(1)+arcTheta(2)),2*pi);
    alph2 = mod(alph1-pi,2*pi);
    tri_proj = [sum((tri-c1).*x1,2),sum((tri-c1).*x2,2)];
    mu = mean(tri_proj);
    N1 = norm([cos(alph1)*R0,sin(alph1)*R0] - mu,2);
    N2 = norm([cos(alph2)*R0,sin(alph2)*R0] - mu,2);
    if(N1<N2)
        alph0 = alph1;
    else
        alph0 = alph2;
    end
    arcTheta = sort(arcTheta);
    if(alph0<arcTheta(1) || alph0>arcTheta(2))
        arcTheta(2) = arcTheta(2)-2*pi;
    end
end
% Add area of circular segments to that of polygon
A = abs(A);
if(~isempty(arcTheta))
    alph = abs(arcTheta(:,1)-arcTheta(:,2));
    A = A+sum(R0^2/2*(alph-sin(alph)));
end

return;

figure(1);
cla
tri_proj = [sum((tri-c1).*x1,2),sum((tri-c1).*x2,2)];
plot(R0*cos(linspace(0,2*pi,1e3)),R0*sin(linspace(0,2*pi,1e3)),'-','color',[0.5,0.5,0.5]);
hold on;
line([tri_proj(1,1),tri_proj(2,1)],[tri_proj(1,2),tri_proj(2,2)],'color',[0.5,0.5,0.5]);
line([tri_proj(3,1),tri_proj(2,1)],[tri_proj(3,2),tri_proj(2,2)],'color',[0.5,0.5,0.5]);
line([tri_proj(1,1),tri_proj(3,1)],[tri_proj(1,2),tri_proj(3,2)],'color',[0.5,0.5,0.5]);
set(gca,'DataAspectRatio',[1,1,1]);

fill(in_proj(:,1),in_proj(:,2),'k','EdgeColor','k','FaceAlpha',0.4,'LineWidth',1);

for i = 1:size(arcTheta,1)
    t = linspace(arcTheta(i,1),arcTheta(i,2),1e2);
    fill(cos(t)*R0,sin(t)*R0,'k','FaceAlpha',0.4,'EdgeColor','k','LineWidth',1);
end
title([num2str(A/(pi*R0^2)*100,3) '%'])
xticks([]); yticks([]);