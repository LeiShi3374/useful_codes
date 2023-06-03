close all; clear all; clc;

fname = 'log_BO_epi.txt';
data = load(fname);
n    = size(data,2);
t    = data(:,1);
V    = data(:,2);

s    = data(:,5);
Ta   = data(:,n);

%% Plots
figure('units','normalized','outerposition',[0.65 0.4 0.6 0.4]);
    %% V, Ta
    subplot(1,2,1);
    plot(t./1000.0, V, 'k-', 'LineWidth', 4);
    hold on;
    yyaxis right;
    plot(t./1000.0, Ta./1000, 'r--', 'LineWidth', 4);
    
    set(gca, 'FontName', 'Times', 'FontSize', 18, 'LineWidth', 2);
    ax = gca;
    ax.XLim = [0.0 1.0];
    xlabel('Time (s)',...
        'FontSize',18,'FontWeight','bold','Color','k');
    
    yyaxis left;
    ax.YLim = [-100 60];
    ylabel('Action Potential (mV)', ...
        'FontSize',18,'FontWeight','bold','Color','k');
    
    yyaxis right;
    ax.YLim = [0 100.0];
    ylabel('Active Stress (KPa)', ...
        'FontSize',18,'FontWeight','bold','Color',[0.8 0 0]);
    yticks([0:20:100]);

    dim = [0.31 0.75 0.25 0.1];
    str = {'Bueno-Orovio'; 'Active Stress'};
    annotation('textbox',dim,'String',str,...
        'FitBoxToText','on','FontSize',20, 'FontWeight','bold',...
        'BackgroundColor','w','EdgeColor','w');
    hold off;

    %% s, Ta
    subplot(1,2,2);
    plot(t./1000.0, s, 'k-', 'LineWidth', 4);
    hold on;
    yyaxis right;
    plot(t./1000.0, Ta./1000, 'r--', 'LineWidth', 4);
    
    set(gca, 'FontName', 'Times', 'FontSize', 18, 'LineWidth', 2);
    ax = gca;
    ax.XLim = [0.0 1.0];
    xlabel('Time (s)',...
        'FontSize',18,'FontWeight','bold','Color','k');
    
    yyaxis left;
    ax.YLim = [0 0.9];
    ylabel('s', ...
        'FontSize',18,'FontWeight','bold','Color','k');
    yticks([0:0.2:1]);
    
    yyaxis right;
    ax.YLim = [0 100.0];
    ylabel('Active Stress (KPa)', ...
        'FontSize',18,'FontWeight','bold','Color',[0.8 0 0]);
    yticks([0:20:100]);

    dim = [0.75 0.75 0.25 0.1];
    str = {'Bueno-Orovio'; 'Active Stress'};
    annotation('textbox',dim,'String',str,...
        'FitBoxToText','on','FontSize',20, 'FontWeight','bold',...
        'BackgroundColor','w','EdgeColor','w');
    hold off;
    
