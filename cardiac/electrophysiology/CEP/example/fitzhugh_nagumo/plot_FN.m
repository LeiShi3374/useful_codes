close all; clear all; clc;

fname = 'log_FN.txt';
data = load(fname);
n  = size(data,2);
t  = data(:,1);
V  = data(:,2);
r  = data(:,3);

figure('units','normalized','outerposition',[0.65 0.4 0.25 0.36]);
    plot(t, V, 'k-', 'LineWidth', 4);
    hold on;
    yyaxis right;
    plot(t, r, 'r--', 'LineWidth', 4);
    
    set(gca, 'FontName', 'Times', 'FontSize', 18, 'LineWidth', 2);
    ax = gca;
    ax.XLim = [0.0 5];
    xlabel('t',...
        'FontSize',24,'FontWeight','bold','Color','k');
    
    yyaxis left;
    ax.YLim = [-0.9 1.2];
    ylabel('$\phi$', 'interpreter', 'latex', ...
        'FontSize',24,'FontWeight','bold','Color','k');
    yticks([-0.5:0.5:1]);
    
    yyaxis right;
    ax.YLim = [-0.5 0.6];
    ylabel('r', ...
        'FontSize',24,'FontWeight','bold','Color',[0.8 0 0]);
    yticks([-0.4:0.2:0.6]);

    dim = [0.57 0.8 0.25 0.1];
    str = {'Fitzhugh-Nagumo'};
    annotation('textbox',dim,'String',str,...
        'FitBoxToText','on','FontSize',20, 'FontWeight','bold',...
        'BackgroundColor','w','EdgeColor','w');
    hold off;
    