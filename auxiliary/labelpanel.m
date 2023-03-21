function A = labelpanel(x,y,str,lowerCase)


if(nargin<4)
	str = upper(str); % CellPress
else
	if(lowerCase)
		str = lower(str); % NPG
	else
		str = upper(str);
	end
end

A = annotation('textbox', [x,y,0.05,0.05],'String',str, 'LineStyle', ...
	'none', 'FontWeight','bold', 'FontSize',8,'Margin',0);