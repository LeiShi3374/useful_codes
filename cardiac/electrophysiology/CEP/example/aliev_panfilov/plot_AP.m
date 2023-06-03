close all; clear all; clc;

fname = 'log_AP.txt';
data = load(fname);
n  = size(data,2);
t  = data(:,1);
V  = data(:,2);
Ta = data(:,n);

figure('units','normalized','outerposition',[0.65 0.4 0.3 0.48]);
    plot(t./1000.0, V, 'k-', 'LineWidth', 4);
    hold on;
    yyaxis right;
    plot(t./1000.0, Ta./1000.0, 'r--', 'LineWidth', 4);
    
    set(gca, 'FontName', 'Times', 'FontSize', 18, 'LineWidth', 2);
    ax = gca;
    ax.XLim = [0.0 1.3];
    xlabel('Time (s)',...
        'FontSize',18,'FontWeight','bold','Color','k');
    
    yyaxis left;
    ax.YLim = [-100 40];
    ylabel('Action Potential (mV)', ...
        'FontSize',18,'FontWeight','bold','Color','k');
    
    yyaxis right;
    ax.YLim = [0 600.0];
    ylabel('Active Stress (KPa)', ...
        'FontSize',18,'FontWeight','bold','Color',[0.8 0 0]);

    dim = [0.65 0.75 0.25 0.1];
    str = {'Aliev-Panfilov'; 'Active Stress'};
    annotation('textbox',dim,'String',str,...
        'FitBoxToText','on','FontSize',20, 'FontWeight','bold',...
        'BackgroundColor','w','EdgeColor','w');
    hold off;
    