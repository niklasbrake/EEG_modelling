x = rand(100,2);
x(1,:) = [0.5,0.5];
d = vecnorm(x-x(1,:),2,2);
f = max(exp(-d*5),1e-9);
f(1) = 1e-9;

k = 10;
I = interp1(cumsum(f)/sum(f),1:length(f),rand(k,1),'next','extrap');
J = rand(100,1)*3<2;
figureNB(5.5,3);
subplot(1,2,1);
    scatter(x(J==0,1),x(J==0,2),3,[0,0,0],'filled','MarkerEdgeColor','k');
    hold on;
    scatter(x(J==1,1),x(J==1,2),3,[0,0,0],'filled','MarkerEdgeColor','k');
    plot(x(1,1),x(1,2),'.r','MarkerSize',10)

    xticks([]);
    yticks([]);
    xlim([0,1]);
    ylim([0,1]);
    axis off;
    axis square;

for i = 1:k
    i0 = I(i);
    line([x(1,1),x(i0,1)],[x(1,2),x(i0,2)],'color','r','LineWidth',0.75);
end

subplot(1,2,2);
for i = 1:k
    i0 = I(i);
    line([x(1,1),x(i0,1)],[x(1,2),x(i0,2)],'color','r','LineWidth',0.75);
end
hold on;
    idcs = intersect(find(J==0),I);
    scatter(x(idcs,1),x(idcs,2),5,[0,0,0],'filled','MarkerEdgeColor','k');
    hold on;
    idcs = intersect(find(J==1),I);
    scatter(x(idcs,1),x(idcs,2),5,[1,1,1],'filled','MarkerEdgeColor','k');
    plot(x(1,1),x(1,2),'.r','MarkerSize',10)

    xticks([]);
    yticks([]);
    xlim([0,1]);
    ylim([0,1]);
    axis off;
    axis square;
% for i = 1:k
%     i0 = I(i);
%     U = x(1,1)-x(i0,1);
%     V = x(1,2)-x(i0,2);
%     Q = quiver(x(i0,1),x(i0,2),U,V,'color','r','LineWidth',0.75);
% end

% x0 = x(1,1)-0*x(I,1);
% y0 = x(1,2)-0*x(I,1);

% U = x(I,1)-x(1,1);
% V = x(I,2)-x(1,2);

% Q = quiver(x0,y0,U,V,'color','r','LineWidth',0.75);