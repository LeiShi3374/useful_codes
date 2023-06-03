    close all;
    clear; clc;

    f1  = 'load.dat';
    f2 = 'interp_out.dat';

    dt = 0.001;
    
    fid = fopen(f1,'r');
    C = textscan(fid,'%f %f','HeaderLines',1);
    fData = cell2mat(C);
    fclose(fid);

    t1 = fData(:,1);
    Q1 = fData(:,2);
    
    Tp = t1(end);
    nTimePerCycle = int32(Tp/dt);
    
    fid = fopen(f2,'r');
    C = textscan(fid,'%f %f %f','HeaderLines',1);
    fData = cell2mat(C);
    fclose(fid);
    
    t2 = fData(:,1);
    Q2 = fData(:,2);
    
    titleFont = 20;
    lineWidth = 2.5;
    axesFont  = 18;
    axesWidth = 1.5;
    
    %% Figure 1
    figure('units','normalized','outerposition',[0.05 0.08 0.25 0.4]);
        plot(t1, Q1, 'k-', 'LineWidth', lineWidth); hold on;
        plot(t2(1:nTimePerCycle), Q2(1:nTimePerCycle), 'r-',...
            'LineWidth', lineWidth); hold off;
        set(gca,'FontSize',axesFont,'LineWidth',axesWidth);

    figure('units','normalized','outerposition',[0.05 0.08 0.25 0.4]);
        plot(t1, Q1, 'k-', 'LineWidth', lineWidth); hold on;
        plot(t2, Q2, 'r-', 'LineWidth', lineWidth); hold off;
        set(gca,'FontSize',axesFont,'LineWidth',axesWidth);

        