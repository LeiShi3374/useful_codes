close all; clear all; clc;

fname = 'log_dcpld.txt';
data = load(fname);
t  = data(:,1);
Ta = data(:,2);

figure('units','normalized','outerposition',[0.65 0.4 0.3 0.48]);
    plot(t./1000.0, Ta./1000.0, 'k-', 'LineWidth', 4);
    
    set(gca, 'FontName', 'Times', 'FontSize', 18, 'LineWidth', 2);
    ax = gca;
    ax.XLim = [0.0 3];
    xlabel('Time (s)',...
        'FontSize',18,'FontWeight','bold','Color','k');
    
    ax.YLim = [0 100];
    ylabel('Active Stress (KPa)', ...
        'FontSize',18,'FontWeight','bold','Color','k');
    yticks([0:20:100]);
    
    dim = [0.54 0.75 0.25 0.1];
    str = {'Decoupled (Pfaller et al.)'; 'Active Stress'};
    annotation('textbox',dim,'String',str,...
        'FitBoxToText','on','FontSize',20, 'FontWeight','bold',...
        'BackgroundColor','w','EdgeColor','w');
    hold off;
    