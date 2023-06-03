close all; clear all; clc;

fname = 'log_TTP_epi.txt';
data = load(fname);
n    = size(data,2);
t    = data(:,1);
V    = data(:,2);

Ca_i = data(:,5);
gf   = data(:,n);

%% Plots
figure('units','normalized','outerposition',[0.65 0.4 0.6 0.4]);
    %% V, Ta
    subplot(1,2,1);
    plot(t./1000.0, V, 'k-', 'LineWidth', 4);
    hold on;
    yyaxis right;
    plot(t./1000.0, gf, 'r-', 'LineWidth', 4);
    
    set(gca, 'FontName', 'Times', 'FontSize', 18, 'LineWidth', 2);
    ax = gca;
    ax.XLim = [0.0 1.0];
    xlabel('Time (s)',...
        'FontSize',18,'FontWeight','bold','Color','k');
    
    yyaxis left;
    ax.YLim = [-100 40];
    ylabel('Action Potential (mV)', ...
        'FontSize',18,'FontWeight','bold','Color','k');
    
    yyaxis right;
    ax.YLim = [-0.095 0.0];
    ylabel('$\gamma_f$', 'Interpreter','latex', ...
        'FontSize',24,'FontWeight','bold','Color',[0.8 0 0]);
    yticks([-0.08:0.02:0.0]);

    dim = [0.27 0.27 0.25 0.1];
    str = {'tenTusscher-Panfilov'; 'Active Strain (RK/RK)'};
    annotation('textbox',dim,'String',str,...
        'FitBoxToText','on','FontSize',20, 'FontWeight','bold',...
        'BackgroundColor','w','EdgeColor','w');
    hold off;

    %% Ca_i, Ta
    subplot(1,2,2);
    plot(t./1000.0, Ca_i.*1000, 'k-', 'LineWidth', 4);
    hold on;
    yyaxis right;
    plot(t./1000.0, gf, 'r-', 'LineWidth', 4);
    
    set(gca, 'FontName', 'Times', 'FontSize', 18, 'LineWidth', 2);
    ax = gca;
    ax.XLim = [0.0 1.0];
    xlabel('Time (s)',...
        'FontSize',18,'FontWeight','bold','Color','k');
    
    yyaxis left;
    ax.YLim = [0 1];
    ylabel('Ca$_i$ (uM)', 'Interpreter','latex', ...
        'FontSize',24,'FontWeight','bold','Color','k');
    yticks([0:0.2:1.0]);
    
    yyaxis right;
    ax.YLim = [-0.095 0.0];
    ylabel('$\gamma_f$', 'Interpreter','latex', ...
        'FontSize',24,'FontWeight','bold','Color',[0.8 0 0]);
    yticks([-0.08:0.02:0.0]);

    dim = [0.71 0.27 0.25 0.1];
    str = {'tenTusscher-Panfilov'; 'Active Strain (RK/RK)'};
    annotation('textbox',dim,'String',str,...
        'FitBoxToText','on','FontSize',20, 'FontWeight','bold',...
        'BackgroundColor','w','EdgeColor','w');
    hold off;
    