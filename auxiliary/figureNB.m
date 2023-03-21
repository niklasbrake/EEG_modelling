function varargout = figureNB(x,y)
    if(nargin==0)
        x=8;
        y=8;
    end
    set(0,'units','centimeters');
    SS = get(0,'screensize');

    fig = figure('color','w','units','centimeters','ToolBar','none');
    fig.Position = [SS(3)/2-x/2,SS(4)/2-y/2,x,y];
    set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[x,y],'Renderer','Painters');
    set(fig, 'DefaultAxesFontName', 'Arial');
    set(fig, 'DefaultTextFontName', 'Arial');

    if(nargout==1)
        varargout{1} = fig;
    else
        varargout = {};
    end