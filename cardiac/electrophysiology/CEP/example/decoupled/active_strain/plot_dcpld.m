close all; clear all; clc;

fname = 'log_dcpld.txt';
data = load(fname);
t  = data(:,1);
gf = data(:,2);

figure('units','normalized','outerposition',[0.65 0.4 0.3 0.48]);
    plot(t./1000.0, gf.*100.0, 'k-', 'LineWidth', 4);
    
    set(gca, 'FontName', 'Times', 'FontSize', 18, 'LineWidth', 2);
    ax = gca;
    ax.XLim = [0.0 3];
    xlabel('Time (s)',...
        'FontSize',18,'FontWeight','bold','Color','k');
    
    ax.YLim = [-16 1];
    ylabel('Fiber shortening (%)', ...
        'FontSize',18,'FontWeight','bold','Color','k');
    yticks([-20:5:0]);
    
    dim = [0.65 0.75 0.25 0.1];
    str = {'Decoupled'; 'Active Strain'};
    annotation('textbox',dim,'String',str,...
        'FitBoxToText','on','FontSize',20, 'FontWeight','bold',...
        'BackgroundColor','w','EdgeColor','w');
    hold off;
    