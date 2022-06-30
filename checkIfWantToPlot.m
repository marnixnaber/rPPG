function [FigH, F] = checkIfWantToPlot(plotYesNo)

if strcmp(plotYesNo, 'No')
    
    figure('visible','off');
    set(gcf, 'PaperUnits', 'inches');
    x_width=15 ;y_width=9.125;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
    
else
    figure();
    set(gcf, 'PaperUnits', 'inches');
    x_width=15 ;y_width=9.125;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
end

end

