folder = 'E:\Research_Projects\004_Propofol\data\resources\cortical_column_Hagen\segmentations';
F = dir(folder); F = F(3:end);


layers = [0,-300,-300,-750,-1000,-1335];
x0(3,:) = [760,layers(4)];
x0(7,:) = [-300,layers(5)];
x0(1,:) = [1150,layers(2)-20];
x0(6,:) = [220,layers(5)];
x0(10,:) = [550,layers(6)+50];
x0(8,:) = [1250,layers(5)];
x0(11,:) = [-130,layers(6)];
x0(9,:) = [950,layers(6)-50];
x0(2,:) = [-100,layers(2)-150];
x0(5,:) = [1075,layers(4)];
x0(4,:) = [-40,layers(4)];

J = [1,6,2,9,3,8,10,5,7,4,11];
for i = 1:11
    x0(i,1) = (J(i)-1)*300;
end

x_min = Inf;
x_max = -Inf;
y_min = Inf;
y_max = -Inf;

fig = figureNB(8,7);
ax = axes('units','centimeters');
L = 6; H = 2.7;
ax.Position = [fig.Position(3)/2-L/2,fig.Position(4)/2-H/2,L,H];
LY = [0; -78; -500; -767; -975; -1250]*1.2;
for i = 1:length(LY)
    line([-1e3,4e3],[1,1]*LY(i),'color',[0.75,0.75,0.75]);
end

for k = 1:size(x0,1)
    morphData = load(fullfile(folder,F(k).name));
    if(isempty(strfind(F(k).name,'I_')))
        clr = clrsPT.qualitative_CM.red;
    else
        clr = clrsPT.qualitative_CM.blue;
    end
    for i = 1:length(morphData.area)
        x = morphData.x(i,:);
        y = morphData.y(i,:);
        z = morphData.z(i,:);
        L = vecnorm(diff([x;y;z]'));
        r = morphData.area(i)/2*pi/L;
        line(x+x0(k,1),z+x0(k,2),'color',clr,'LineWidth',1./(1+7*exp(-(r-1)/2)));
        x_min = min(x_min,min(x+x0(k,1)));
        x_max = max(x_max,max(x+x0(k,1)));
        y_min = min(y_min,min(z+x0(k,2)));
        y_max = max(y_max,max(z+x0(k,2)));
    end
    set(gca,'DataAspectRatio',[1,1,1])
    drawnow;
end
xlim([x_min,x_max]);
ylim([y_min,y_max]);
axis off;

ty = diff(LY)/2+LY(1:end-1);
str = {'L1','L2/3','L4','L5','L6'};
for i =1 :length(ty)
    text(x_min-(x_max-x_min)/L*0.04,ty(i),str{i},'FontSize',6,'color',[1,1,1]*0.6,'HorizontalAlignment','right')
end



neuronTypes = {'L23E_oi24rpy1';'L23I_oi38lbc1';'L23I_oi38lbc1';'L4E_53rpy1';'L4E_j7_L4stellate';'L4E_j7_L4stellate';'L4I_oi26rbc1';'L4I_oi26rbc1';'L5E_oi15rpy4';'L5E_j4a';'L5I_oi15rbc1';'L5I_oi15rbc1';'L6E_51_2a_CNG';'L6E_oi15rpy4';'L6I_oi15rbc1';'L6I_oi15rbc1'};
nrnAbundance = [26.8,3.2,4.3,9.5,9.5,9.5,5.6,1.5,4.9,1.3,0.6,0.8,14,4.6,1.9,1.9]/100;
mTypes = unique(neuronTypes);
mTypeCount = accumarray(findgroups(neuronTypes),nrnAbundance)

for i = 1:length(mTypes)
    fprintf('%s: %.1f\n',mTypes{i},1e2*mTypeCount(i));
end