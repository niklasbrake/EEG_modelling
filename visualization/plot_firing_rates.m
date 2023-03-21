% figureNB;
mTypes = unique(morphs);
for i = 1:11
    if(isempty(strfind(mTypes{i},'I_')))
        clrs(i,:) = [1,0,0];
    else
        clrs(i,:) = [0,0,1];
    end
    % boxplotNB(i,NN{i},clr,10);
end

m2 = cellfun(@(x)strrep(x,'_','\_'),unique(morphs),'UniformOutput',false);


figureNB(5.8,5.22);
    violin(NN,1:11,clrs,true);
    view([90,90])
    xticks(1:11)
    xticklabels(m2)
    ylim([-1,12])
    ylabel('Firing Frequency (Hz)')
gcaformat(gcf);