% This script produces Figure 1 where the proposed implementation of GCSD is 
% compared to the conventional one


clear
clc
close all

load(fullfile("simulated","sim6.mat"))
data = Data{11};
fs = 0.5;

gcpsd_1 = gcpsd(data,size(data,1));
gcpsd_4 = gcpsd(data,size(data,1)/4);
gcpsd_comb = gcpsd_new(data);

f = linspace(0,fs/2,length(gcpsd_1));
pos = 53;
gcpsd_1 = gcpsd_1(1:pos);
gcpsd_4 = gcpsd_4(1:pos);
gcpsd_comb = gcpsd_comb(1:pos);
f = f(1:pos);

h = tight_subplot(3,1,[0.05,0.05],[0.13,0.05],[0.11,0.05]);
axes(h(1))
plot(f,gcpsd_1,'LineWidth',1.5)
xticklabels([])
ax = gca;
ax.FontSize = 10; 
ylabel('W = N','FontSize',14,'Position',[-0.008,0.05])
xlim([0,0.1])
axes(h(2))
plot(f,gcpsd_4,'LineWidth',1.5)
xticklabels([])
ax = gca;
ax.FontSize = 10; 
ylabel('W = N/4','FontSize',14,'Position',[-0.008,0.0075])
xlim([0,0.1])
axes(h(3))
plot(f,gcpsd_comb,'LineWidth',1.5)
xticks([0,0.015,0.03,0.05,0.065,0.1])
ax = gca;
ax.FontSize = 10; 
ylabel('Proposed','FontSize',14,'Position',[-0.008,0.2])
xlabel('Frequency (Hz)','FontSize',14)
xlim([0,0.1])


exportgraphics(gcf,'C:\Users\Hp\Desktop\Khalifa\PhD\Dissertation\MSWD_paper_new\paper_plots\gcpsd.png','Resolution',1000)