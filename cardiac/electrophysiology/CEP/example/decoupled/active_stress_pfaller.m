clc;
clear;
close all;

% Ventricle
sigma_0 = 1e5;
t_sys   = 0.17;
t_dia   = 0.484;

% Atrium
% sigma_0 = 972.0;
% t_sys   = 0.07;
% t_dia   = 0.14;

gamma     = 0.005;
alpha_min = -30;
alpha_max = 5;

h  = 0.001;
t  = 0:h:1;
nt = length(t);

tau = zeros(2,nt);

%% Euler
for i = 1 : length(t)-1
    Sp = 1/2*(1+tanh((t(i)-t_sys)/gamma));
    Sn = 1/2*(1-tanh((t(i)-t_dia)/gamma));
    f = Sp.*Sn;
    a = alpha_max*f + alpha_min*(1-f);
    k1 = -abs(a)*tau(1,i) + sigma_0*max(0,a);
    tau(1,i+1) = tau(1,i) + h*k1;
end

%% Runge-Kuta
for i = 1 : length(t)-1
    Sp = 1/2*(1+tanh((t(i)-t_sys)/gamma));
    Sn = 1/2*(1-tanh((t(i)-t_dia)/gamma));
    f = Sp.*Sn;
    a = alpha_max*f + alpha_min*(1-f);
    k1 = -abs(a)*tau(2,i) + sigma_0*max(0,a);

    Sp = 1/2*(1+tanh((t(i)+0.5*h-t_sys)/gamma));
    Sn = 1/2*(1-tanh((t(i)+0.5*h-t_dia)/gamma));
    f = Sp.*Sn;
    a = alpha_max*f + alpha_min*(1-f);
    k2 = -abs(a)*(tau(2,i) + 0.5*h*k1) + sigma_0*max(0,a);

    Sp = 1/2*(1+tanh((t(i)+0.5*h-t_sys)/gamma));
    Sn = 1/2*(1-tanh((t(i)+0.5*h-t_dia)/gamma));
    f = Sp.*Sn;
    a = alpha_max*f + alpha_min*(1-f);
    k3 = -abs(a)*(tau(2,i) + 0.5*h*k2) + sigma_0*max(0,a);

    Sp = 1/2*(1+tanh((t(i)+h-t_sys)/gamma));
    Sn = 1/2*(1-tanh((t(i)+h-t_dia)/gamma));
    f = Sp.*Sn;
    a = alpha_max*f + alpha_min*(1-f);
    k4 = -abs(a)*(tau(2,i) + h*k3) + sigma_0*max(0,a);
    tau(2,i+1) = tau(2,i) + 1/6*h*(k1 + 2*k2 + 2*k3 + k4);
end

%% Plot
lineWidth = 2.5;
axesFont  = 18;
axesWidth = 1.5;

figure('units','normalized','outerposition',[0.05 0.08 0.25 0.4]);
plot(t, tau(1,:)./1000.0,'k', 'LineWidth',4)
hold on
plot(t, tau(2,:)./1000.0,'r--', 'LineWidth',4)
hold off

set(gca,'FontSize',axesFont,'LineWidth',axesWidth);
xlabel('time','FontSize',axesFont);
ylabel('Active Stress (KPa)','FontSize',axesFont);
legend('Euler','RK4','FontSize',axesFont,'Location','northeast');
legend('boxoff');

%lb = linspace(1,length(t),501);
%plot(t(lb), tau(lb),'k','LineWidth',4)


%writematrix([length(lb), 256; [t(lb)', tau(lb)'*10]],'stress.dat','Delimiter',',')

