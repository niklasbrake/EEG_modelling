function gcaformat(graphicsObject,article)
    if(nargin<2)
        article = true;
    end
    if(nargin==0)
        nbFormat(gca,article)
    elseif(strcmp(graphicsObject.Type,'axes'))
        nbFormat(graphicsObject,article)
    elseif(strcmp(graphicsObject.Type,'figure'))
        axs = get(graphicsObject,'children');
        for i = 1:length(axs)
            if(strcmp(axs(i).Type,'axes'))
                nbFormat(axs(i),article);
            end
        end
    else
        error('invalid input type');
    end

end
function nbFormat(ax,article)
    if(article)
        nbFormat_article(ax);
    else
        nbFormat_presentation(ax);
    end
end
function nbFormat_article(ax)
    set(ax,'box','off');
    set(ax,'tickdir','out');
    % set(ax,'linewidth',0.75);
    set(ax,'linewidth',0.5);
    set(ax,'fontsize',7);
    set(ax.Title,'FontSize',7);
    xax = get(ax,'xaxis');
    xax.Label.FontSize = 7;
    xax.TickLabelRotation = 0;
    yax = get(ax,'yaxis');
    for i = 1:length(yax)
        yax(i).Label.FontSize = 7;
        yax(i).TickLabelRotation = 0;
    end
    tickLength = 0.05; % 1/2 mm
    U = ax.Units;
    set(ax,'Units','centimeters');
    L = max(ax.Position(3:4));
    t = tickLength/L;
    set(ax,'TickLength',[t,t]);
    set(ax,'Units',U);
end
function nbFormat_presentation(ax)
    set(ax,'box','off');
    set(ax,'tickdir','out');
    set(ax,'linewidth',1.5);
    set(ax,'fontsize',18);
    set(ax.Title,'FontSize',22);
    xax = get(ax,'xaxis');
    xax.Label.FontSize = 22;
    yax = get(ax,'yaxis');
    for i = 1:length(yax)
        yax(i).Label.FontSize = 22;
    end
    tickLength = 0.2; % 2 mm
    U = ax.Units;
    set(ax,'Units','centimeters');
    L = max(ax.Position(3:4));
    t = tickLength/L;
    set(ax,'TickLength',[t,t]);
    set(ax,'Units',U);
end