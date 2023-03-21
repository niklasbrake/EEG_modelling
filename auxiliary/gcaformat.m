function gcaformat(graphicsObject)
    if(nargin==0)
        nbFormat(gca)
    elseif(strcmp(graphicsObject.Type,'axes'))
        nbFormat(graphicsObject)
    elseif(strcmp(graphicsObject.Type,'figure'))
        axs = get(graphicsObject,'children');
        for i = 1:length(axs)
            if(strcmp(axs(i).Type,'axes'))
                nbFormat(axs(i));
            end
        end
    else
        error('invalid input type');
    end

end
function nbFormat(ax)
    set(ax,'box','off');
    set(ax,'tickdir','out');
    set(ax,'linewidth',0.75);
    set(ax,'fontsize',7);
    set(ax.Title,'FontSize',7);
    xax = get(ax,'xaxis');
    xax.Label.FontSize = 7;
    yax = get(ax,'yaxis');
    for i = 1:length(yax)
        yax(i).Label.FontSize = 7;
    end
    tickLength = 0.05; % 1/2 mm
    U = ax.Units;
    set(ax,'Units','centimeters');
    L = max(ax.Position(3:4));
    t = tickLength/L;
    set(ax,'TickLength',[t,t]);
    set(ax,'Units',U);
end