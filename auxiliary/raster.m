function varargout = raster(trialNumber,spikeTime,fig)
    if(nargin<2)
        spikeTime = 0;
    end
    if(nargin<3)
        figureNB;
    else
        figure(fig);
    end
    if(sum(size(spikeTime)>1)>1)
        A = spikeTime;
        numspikes = sum(spikeTime(:));
        spikeTime = zeros(numspikes,1);
        trialNumber = zeros(size(spikeTime));
        j = 0;
        for i = 1:size(A,1)
            temp = find(A(i,:));
            m = length(temp);
            spikeTime(j+1:j+m) = temp;
            trialNumber(j+1:j+m) = i+0*temp;
            j = j+m;
        end
    end
    numspikes = length(trialNumber);
    xx=nan(3*numspikes,1);
    yy=nan(3*numspikes,1);
    yy(1:3:3*numspikes)=trialNumber;
    yy(2:3:3*numspikes)=trialNumber+15;
    xx(1:3:3*numspikes)=spikeTime;
    xx(2:3:3*numspikes)=spikeTime;

    h = plot(xx,yy,'-k','LineWidth',1);
    xlim([min(spikeTime),max(spikeTime)])
    ylim([min(trialNumber),max(trialNumber)+1])
    gcaformat;
    xlabel('Time');
    set(get(gca,'yaxis'),'visible','off')


    if(nargout==1)
        varargout{1} = h;
    else
        varargout = {};
    end