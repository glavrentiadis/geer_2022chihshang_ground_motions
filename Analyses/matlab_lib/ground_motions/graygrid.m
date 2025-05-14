function [] = gridgray ()

% work around because cannot set gridlines color and tick label color separatly

% set grid lines to lines not dotted
set (gca,'minorgridlinestyle','-')
set (gca,'gridlinestyle','-')

grid on;

% set gridlines to gray
set (gca,'xcolor',[.7 .7 .7])
set (gca,'ycolor',[.7 .7 .7])

% now write tick labels in black over the light gray
c_axes = copyobj (gca,gcf); 
set (c_axes, 'color', 'none', 'xcolor', 'k', 'xgrid', 'off', 'ycolor','k', 'ygrid','off');












