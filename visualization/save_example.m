folder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\test\triangle_int';
i = 1;
while(exist(fullfile(folder,[int2str(i) '.mat'])))
    i = i+1;
end
save(fullfile(folder,[int2str(i) '.mat']),'tri');