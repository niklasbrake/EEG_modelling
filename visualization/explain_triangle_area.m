c0 = [0,0,0];

R = 1;

% tri = [0.5,-1.5,-0.5; -1.1,-0.5,-0.6; -0.5,0,0.5]+[0,1,0];
tri = randn(3);
% folder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\test\triangle_int';
% load(fullfile(folder,'8.mat'));
params = null([tri,ones(3,1)])';
N = params(1:3);
rho = dot(N,c0)-params(4)/norm(N,2);
% Case 1
if(abs(rho)>R)
    A = 0
    return;
end

pairs = [1,2;2,3;3,1];
A = sum(tri(pairs(:,1),:).^2,2) - R^2;
C = sum((tri(pairs(:,1),:)-tri(pairs(:,2),:)).^2,2);
B = sum(tri(pairs(:,2),:).^2,2) - A - C - R^2;
in = [];
T = [];
for k = 1:size(pairs,1)
    t = roots([C(k),B(k),A(k)]);
    if(sum(t<0) == 2 || sum(t>1) == 2)
        t = t*sqrt(-1);
    end
    T = [T;pairs(k,:),t(1);pairs(k,:),t(2)];
    in = [in;tri(pairs(k,1),:).*(1-t) + tri(pairs(k,2),:).*t];
end

in(T(:,3)<0,:) = tri(T(T(:,3)<0,1),:);
in(T(:,3)>1,:) = tri(T(T(:,3)>1,2),:);
in(imag(T(:,3))~=0,:) = [];
T(imag(T(:,3))~=0,:) = [];
[T,idcs] = sortrows(T);
in = in(idcs,:);

R0 = sqrt(R^2-rho^2);
c1 = rho*N/norm(N,2);
[x1,x2] = getOrthBasis(params(1:3));
x2 = x2;
x1 = x1;
tri_proj = [sum((tri-c1).*x1,2),sum((tri-c1).*x2,2)];
in_proj = [sum((in-c1).*x1,2),sum((in-c1).*x2,2)];

fig = figureNB;
fig.Position = [14.7955 7.3025 6.7204 12.9117];
subplot(1,2,1);
    [xs,ys,zs] = sphere(30);
    ms = mesh(xs*R,ys*R,zs*R); hold on;
    ms.LineStyle = ':';
    ms.EdgeColor = 1-[1,1,1]*0.6;
    ms.FaceColor = 1-[1,1,1]*0.1;
    ms.FaceAlpha = 0.8;
    set(gca,'DataAspectRatio',[1,1,1]);
    patch([tri(1,1);tri(2,1);tri(3,1)],[tri(1,2);tri(2,2);tri(3,2)],[tri(1,3);tri(2,3);tri(3,3)],'b','FaceAlpha',0.5);
    hold on;
    for i = 1:3
        plot3(tri(i,1),tri(i,2),tri(i,3),'.k','MarkerSize',20)
        text(tri(i,1)-0.3,tri(i,2),tri(i,3),int2str(i))
    end
    for i = 1:size(in,1)
        plot3(in(i,1),in(i,2),in(i,3),'.r','MarkerSize',10);
    end
    lighting GOURAUD
    % plot3(0,0,0,'.k','MarkerSize',20);
    % line([0,c1(1)],[0,c1(2)],[0,c1(3)],'color','r');
subplot(1,2,2);
    patch([tri_proj(1,1);tri_proj(2,1);tri_proj(3,1)],[tri_proj(1,2);tri_proj(2,2);tri_proj(3,2)],'b','FaceAlpha',0.5);
    hold on;
    set(gca,'DataAspectRatio',[1,1,1]);
    theta = linspace(0,2*pi,1e3);
    fill(R0*cos(theta),R0*sin(theta),'w');
    xl = get(gca,'xlim');
    yl = get(gca,'ylim');
    xlim(xl); ylim(yl)
    axis off;
    % Cases 2
    if(isempty(in))
        if(norm(mean(tri_proj),2)<R0)
            A = pi*R0^2
        else
            A = 0
        end
        title([num2str(A/(pi*R0^2)*100,3) '%'])
        return
    end
    for i = 1:3
        plot(tri_proj(i,1),tri_proj(i,2),'.k','MarkerSize',20); hold on;
        text(tri_proj(i,1)-0.3,tri_proj(i,2),int2str(i))
    end

    % Cases 3-9
    % for i = 1:size(in_proj,1)
    %     plot(in_proj(i,1),in_proj(i,2),'.r','MarkerSize',10);
    % end
    [theta,r] = cart2pol(in_proj(:,1),in_proj(:,2));
    theta = mod(theta,2*pi);
    theta = [theta;theta(1)];
    r = [r;r(1)];



    x = cos(theta).*r;
    y = sin(theta).*r;
    A = sum((x(2:end)+x(1:end-1)).*(diff(y)))/2;
    clockwise = (A<0);
    A = abs(A);
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
    % Case 3
    if(size(T,1)==2)
        arcTheta = mod(arcTheta,2*pi);
        alph1 = mod(0.5*(arcTheta(:,1)+arcTheta(:,2)),2*pi);
        alph2 = mod(alph1-pi,2*pi);
        mu = mean(tri_proj);
        N1 = norm([cos(alph1)*R0,sin(alph1)*R0] - mu,2);
        N2 = norm([cos(alph2)*R0,sin(alph2)*R0] - mu,2);
        if(N1<N2)
            alph0 = alph1;
        else
            alph0 = alph2;
        end
        arcTheta = sort(arcTheta);
        if(alph0>arcTheta(1) & alph0<arcTheta(2))
            % do nothing
        else
            arcTheta(2) = arcTheta(2)-2*pi;
        end
        % return;
    end
    if(~isempty(arcTheta))
        alph = abs(arcTheta(:,1)-arcTheta(:,2));
        A = A+sum(R0^2/2*(alph-sin(alph)));
    end

    fill(x,y,'r','EdgeColor','k','FaceAlpha',0.4,'LineWidth',1);
    for i = 1:size(arcTheta,1)
        t = linspace(arcTheta(i,1),arcTheta(i,2),1e2);
        fill(cos(t)*R0,sin(t)*R0,'r','FaceAlpha',0.4,'EdgeColor','k','LineWidth',1);
    end
    % title([num2str(A/(pi*R0^2)*100,3) '%'])
    xticks([]);
    yticks([]);
