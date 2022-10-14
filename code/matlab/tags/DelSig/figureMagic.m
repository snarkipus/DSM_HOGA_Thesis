function figureMagic( xRange, dx, xLab, yRange, dy, yLab, size, name )
%figureMagic( xRange, dx, xLab, yRange, dy, yLab, size, name )
if nargin >=8
    set(gcf,'MenuBar','none'); 
    set(gcf,'NumberTitle','off'); 
    set(gcf,'Name',name);
end
if nargin >=7 & ~isempty(size)
    set(gcf,'PaperUnits','inches','PaperPosition', [0.5 0.5 0.5+size]);
end
if nargin <6	% Defaulted yRange etc
    v = axis; 
    yRange = v(3:4);
    DoYaxis = 0;
else
    DoYaxis = 1;
end
axis([xRange yRange] );
grid on; xlabel(''); ylabel('');
set(gca,'FontName','helvetica');
set(gca,'FontSize', 11);
set(gca, 'Box', 'on');
if ~isempty(dx)
    xtix = xRange(1):dx:xRange(2);
    set(gca,'XTick', xtix);
    set(gca,'XTickLabel', axisLabels(xtix,xLab));
end
if DoYaxis
    ytix = yRange(1):dy:yRange(2);
    set(gca,'YTick', ytix);
    set(gca,'YTickLabel', axisLabels(ytix,yLab));
end
