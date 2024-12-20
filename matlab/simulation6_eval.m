clear
clc
close all


%% Sim 6
save_plots = 1;
casee = 6;
TR = 2;                                                   % Repetition time 
fs = 1/TR;                                                % Sampling frequency
t = 0:1/fs:668-1/fs;
f1 = 0.05;
f2 = 0.03;
N = 100;                                                 % number of repetition (realizations)
freq_ranges = [0.01,0.015,0.03,0.05,0.065,0.1];
min_freq = 0;                                        % window sizes for the Windowed Phase Sync. Measures

true_plv = zeros(1,length(freq_ranges));
delphi1 = 2*pi./(1+exp(-0.01*(t-334)));
delphi2 = pi/4;

x_plv1 = cos(2*pi*f1*t);
y_plv1 = cos(2*pi*f1*t + delphi1);
x_plv2 = cos(2*pi*f2*t);
y_plv2 = cos(2*pi*f2*t + delphi2);
phase_x1 = angle(hilbert(x_plv1)); phase_y1 = angle(hilbert(y_plv1));
dPhi1 = phase_x1 - phase_y1;
e1 = exp(1i*dPhi1);
phase_x2 = angle(hilbert(x_plv2)); phase_y2 = angle(hilbert(y_plv2));
dPhi2 = phase_x2 - phase_y2;
e2 = exp(1i*dPhi2);
pos1 = find(freq_ranges == f1);
pos2 = find(freq_ranges == f2);
true_plv(pos1) = abs(mean(e1));
true_plv(pos2) = abs(mean(e2));
true_cosdelphi1 = cos(phase_x1 - phase_y1);
true_cosdelphi = true_cosdelphi1;
% filename = fullfile(pwd,'plots',['sim',num2str(casee) , '_ground_truth.png']);
% simulation2_groundTruthPlot(true_cosdelphi,t,filename)

load(fullfile("decomposed","sim6_MVMD.mat"))
load(fullfile("decomposed","sim6_MSWD.mat"))


%% MVMD
for m = 1:N
   imf = imfs_MVMD{m};
    if m == 1
        plv_length = length(freq_ranges);
    else
        plv_length = length(plv(m-1,:));
    end
    [COSDELPHI11,plv1,mfreq1] = phase_sync_analysis12356(imf,...
    "MVMD",fs,min_freq,plv_length,f1);
    COSDELPHI1(m,:) = COSDELPHI11;
    mfreq{m} = mfreq1;
    plv(m,1:length(plv1)) = plv1;
end
[plv_MVMD,counts_MVMD] = sort_values(mfreq,freq_ranges,plv);
[mean_plv_MVMD,lower_plv_MVMD,upper_plv_MVMD] = confidence_interval(plv_MVMD);

cosdelphi_mvmd = COSDELPHI1;

disp('MVMD')
clear COSDELPHI1 CCORSW mfreq plv


%% MSWD
for m = 1:N
    imf = imfs_MSWD{m};
    if m == 1
        plv_length = length(freq_ranges);
    else
        plv_length = length(plv(m-1,:));
    end
    [COSDELPHI11,plv1,mfreq1] = phase_sync_analysis12356(imf,...
    "MSWD",fs,min_freq,plv_length,f1);
    COSDELPHI1(m,:) = COSDELPHI11;
    mfreq{m} = mfreq1;
    plv(m,1:length(plv1)) = plv1;
end

[plv_MSWD,counts_MSWD] = sort_values(mfreq,freq_ranges,plv);
[mean_plv_MSWD,lower_plv_MSWD,upper_plv_MSWD] = confidence_interval(plv_MSWD);

cosdelphi_mswd = COSDELPHI1;


%% PLOTS
filename = fullfile(pwd,'plots',['sim', num2str(casee), '_TVPS.png']);
if save_plots
    simulation12356Plots(true_cosdelphi,cosdelphi_mvmd,cosdelphi_mswd,t,filename)
else
    simulation12356Plots(true_cosdelphi,cosdelphi_mvmd,cosdelphi_mswd,t)
end
disp('MSWD')


E_plv_MVMD = mean_plv_MVMD - lower_plv_MVMD;
E_plv_MSWD = mean_plv_MSWD - lower_plv_MSWD;

plv_min = min([min(lower_plv_MVMD),...
    min(lower_plv_MSWD),min(true_plv)]);
plv_max = max([max(upper_plv_MVMD),...
    max(upper_plv_MSWD),max(true_plv)]);


figure
errorbar(1:size(plv_MVMD,2),mean_plv_MVMD,E_plv_MVMD,'bo',"LineStyle","none",'LineWidth',1)
hold on
errorbar(size(plv_MVMD,2)+1:size(plv_MVMD,2)+size(plv_MSWD,2),mean_plv_MSWD,E_plv_MSWD,'go',"LineStyle","none",'LineWidth',1)
plot([true_plv,true_plv],'ro','LineWidth',1)
% title('Simulation 7: PLV')
% xlabel('Method')
ylabel('PLV','FontSize',16)
ylim([0,plv_max+0.02])
xlim([1,size(plv_MVMD,2)+size(plv_MSWD,2)])
start = (1 + length(freq_ranges))/2;
stop = 2*length(freq_ranges) + 1 - start;
ticks = start:length(freq_ranges):stop;
xticks(ticks)
xticklabels({'MVMD-PS','MSwD-PS'})
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',16)
xline((length(freq_ranges):length(freq_ranges):length(freq_ranges))+0.5,'--')
hold off
if save_plots
    exportgraphics(gcf,fullfile(pwd,'plots',['sim',num2str(casee),'_plv.png']),'Resolution',1000)
end