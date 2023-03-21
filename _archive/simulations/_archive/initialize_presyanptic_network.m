function initialize_presyanptic_network(spikingFile,tmax,N,m);
    [neuronID,spikeTime,ei,B] = simulatespikes(tmax,N,m);
    fid = fopen(spikingFile,'w');
    [uM,~,nTemp] = unique(neuronID(:));
    I = accumarray(nTemp,1:numel(nTemp),[],@(x){sort(x)});
    UC = 0;
    flop = rand(max(neuronID),1);
    for i = 1:max(neuronID)
        if(uM(i-UC)==i)
            ts = spikeTime(I{i-UC});
        else
            ts = [];
            UC = UC+1;
        end
        if(ei(i)==0)
            if(flop(i)<0.5)
                target = 'eA';
            elseif(flop(i)<=1)
                target = 'eB';
            end
        else
            if(flop(i) < 0.5)
                target = 'iA';
            elseif(flop(i) <= 0.9)
                target = 'iB';
            elseif(flop(i) <= 1)
                target = 'iS';
            end
        end
        s = sprintf('%0.1f,',sort(ts));
        fprintf(fid,'%s,%s\n',target,s(1:end-1));
    end
    fclose(fid);