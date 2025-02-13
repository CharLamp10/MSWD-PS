clear
clc
close all

TR = 2;
fs = 1/TR;
t = 0:1/fs:668-1/fs;
f = 0.05;
N = 100;                                         
SNR = -3;


delphi1 = 4*pi/334.*(t-334).*(t-334>=0);
delphi2 = 2*pi./(1+exp(-0.01*(t-334)));

x1 = cos(2*pi*f*t);
y1 = cos(2*pi*f*t + delphi1);
phase_x1 = angle(hilbert(x1)); phase_y1 = angle(hilbert(y1));
true_cosdelphi1 = cos(phase_x1 - phase_y1);

x2 = cos(2*pi*f*t);
y2 = cos(2*pi*f*t + delphi2);
phase_x2 = angle(hilbert(x2)); phase_y2 = angle(hilbert(y2));
true_cosdelphi2 = cos(phase_x2 - phase_y2);

figure
h = tight_subplot(2, 1, [0.15 0.02],[0.12 0.1],[0.12 0.05]);
axes(h(1))
plot(t,true_cosdelphi1,"Color",'r','LineWidth',2)
title('Simulation 1 Ground Truth CRP', 'FontSize', 18)
ax = gca;
ax.FontSize = 11; 
ylabel('CRP','FontSize',16)
xlim([0,max(t)])
ylim([-1,1])
grid on

axes(h(2))
plot(t,true_cosdelphi2,"Color",'r','LineWidth',2)
title('Simulations 2,3,5,6 Ground Truth CRP', 'FontSize', 18)
ax = gca;
ax.FontSize = 11; 
ylabel('CRP','FontSize',16)
xlabel('Time (sec)','FontSize',16)
xlim([0,max(t)])
ylim([-1,1])
grid on

saveas(gcf,fullfile('plots_new','simulations12356_GT.png'))