clear
clc
close all

%% Sim 7
save_plots = 1;
casee = 7;
TR = 2;         % Repetition time
fs = 1/TR;      % Sampling frequency
t = 0:1/fs:668-1/fs;
f = 0.05;
N = 100;        % number of repetition (realizations)
freq_ranges = 0.01:0.02:0.1;
min_freq = 0;                                       

true_plv = zeros(1,length(freq_ranges));
x1 = cos(2*pi*f*t);
y1 = cos(2*pi*f*t - pi/4);
z1 = cos(2*pi*f*t - pi/2);
x2 = cos(2*pi*f*t + pi);
y2 = cos(2*pi*f*t + 3*pi/4);
z2 = cos(2*pi*f*t + pi/2);

cosdelphi_grd(:,1) = repmat(cos(pi/4),[1,length(t)]); cosdelphi_grd(:,2) = repmat(cos(pi/2),[1,length(t)]);
cosdelphi_grd(:,3) = repmat(cos(-pi),[1,length(t)]); cosdelphi_grd(:,4) = repmat(cos(-3*pi/4),[1,length(t)]);
cosdelphi_grd(:,5) = repmat(cos(-pi/2),[1,length(t)]); cosdelphi_grd(:,6) = repmat(cos(-pi/4 + pi/2),[1,length(t)]);
cosdelphi_grd(:,7) = repmat(cos(-pi/4-pi),[1,length(t)]); cosdelphi_grd(:,8) = repmat(cos(-pi/4 - 3*pi/4),[1,length(t)]);
cosdelphi_grd(:,9) = repmat(cos(-pi/4-pi/2),[1,length(t)]); cosdelphi_grd(:,10) = repmat(cos(-pi/2 - pi),[1,length(t)]);
cosdelphi_grd(:,11) = repmat(cos(-pi/2-3*pi/4),[1,length(t)]); cosdelphi_grd(:,12) = repmat(cos(-pi/2- pi/2),[1,length(t)]);
cosdelphi_grd(:,13) = repmat(cos(pi - 3*pi/4),[1,length(t)]); cosdelphi_grd(:,14) = repmat(cos(pi-pi/2),[1,length(t)]);
cosdelphi_grd(:,15) = repmat(cos(3*pi/4 - pi/2),[1,length(t)]); 
pwelch_window = length(x1);

phase_x1 = angle(hilbert(x1)); 
phase_y1 = angle(hilbert(y1));
phase_z1 = angle(hilbert(z1));

phase_x2 = angle(hilbert(x2)); 
phase_y2 = angle(hilbert(y2));
phase_z2 = angle(hilbert(z2));

pos = find(freq_ranges == f);
true_plv(pos) = 1;

indx = nchoosek(1:6,2);

load(fullfile("decomposed","sim7_MVMD.mat"))
load(fullfile("decomposed","sim7_MSWD.mat"))


%% MVMD
for m = 1:N
    imf = imfs_MVMD{m};
    imf = squeeze(imf);
    if m == 1
        plv_length = [length(freq_ranges),size(indx,1)];
    else
        plv_length = [length(plv(m-1,:,1)),size(indx,1)];
    end
    [COSDELPHI,plv1,mfreq1] = phase_sync_analysis47(imf,"MVMD",indx,...
     fs,f,min_freq,plv_length);
    COSDELPHI_MVMD{m} = COSDELPHI;
    for i = 1:size(COSDELPHI,2)
        nrmse_MVMD(m,i) = rmse(COSDELPHI(:,i),cosdelphi_grd(:,i));
    end
    plv(m,1:size(plv1,1),1:size(plv1,2)) = plv1;
    mfreq{m} = mfreq1;
end

for i = 1:size(indx,1)
    [plv_MVMD(:,:,i),counts_MVMD(:,:,i)] = ...
        sort_values(mfreq,freq_ranges,plv(:,:,i));
end  
plv_MVMD = mean(plv_MVMD,3);
counts_MVMD = mean(counts_MVMD,3);
[mean_plv_MVMD,lower_plv_MVMD,upper_plv_MVMD] = confidence_interval(plv_MVMD);

filename1 = fullfile(pwd,'plots',['sim', num2str(casee), '_MVMD']);
if save_plots
    simulation7Plots(COSDELPHI_MVMD,cosdelphi_grd,[],[],t,filename1);
else
    simulation7Plots(COSDELPHI_MVMD,cosdelphi_grd,[],[],t);
end
disp('MVMD')
clear mfreq plv


%% MSWD
for m = 1:N
    imf = imfs_MSWD_CL{m};
    if m == 1
        plv_length = [length(freq_ranges),size(indx,1)];
    else
        plv_length = [length(plv(m-1,:,1)),size(indx,1)];
    end
    [COSDELPHI,plv1,mfreq1] = phase_sync_analysis47(imf,"MSWD",indx,...
     fs,f,min_freq,plv_length);
    COSDELPHI_MSWD{m} = COSDELPHI;
    for i = 1:size(COSDELPHI,2)
        nrmse_MSWD(m,i) = rmse(COSDELPHI(:,i),cosdelphi_grd(:,i));
    end
    plv(m,1:size(plv1,1),1:size(plv1,2)) = plv1;
    mfreq{m} = mfreq1;
end
for i = 1:size(indx,1)
    [plv_MSWD(:,:,i),counts_MSWD(:,:,i)] = ...
        sort_values(mfreq,freq_ranges,plv(:,:,i));
end
plv_MSWD = mean(plv_MSWD,3);
counts_MSWD = mean(counts_MSWD,3);
[mean_plv_MSWD,lower_plv_MSWD,upper_plv_MSWD] = confidence_interval(plv_MSWD);


%% MSWD PLOTS
filename1 = fullfile(pwd,'plots',['sim', num2str(casee), '_MSWD']);
filename2 = fullfile(pwd,'plots',['sim', num2str(casee), '_RMSE']);
if save_plots
    simulation7Plots(COSDELPHI_MSWD,cosdelphi_grd,mean(nrmse_MVMD),mean(nrmse_MSWD),t,filename1,filename2);
else
    simulation7Plots(COSDELPHI_MSWD,cosdelphi_grd,mean(nrmse_MVMD),mean(nrmse_MSWD),t);
end
disp('MSWD')


E_plv_MVMD = mean_plv_MVMD - lower_plv_MVMD;
E_plv_MSWD = mean_plv_MSWD - lower_plv_MSWD;

plv_min = min([min(lower_plv_MVMD),...
    min(lower_plv_MSWD),min(true_plv)]);
plv_max = max([max(upper_plv_MVMD),...
    max(upper_plv_MSWD),max(true_plv)]);

figure
errorbar(1:size(plv_MVMD,2),mean_plv_MVMD,E_plv_MVMD,'bo',"LineStyle","none")
hold on
errorbar(size(plv_MVMD,2)+1:size(plv_MVMD,2)+size(plv_MSWD,2),mean_plv_MSWD,E_plv_MSWD,'go',"LineStyle","none")
plot([true_plv,true_plv],'ro')
% title('Simulation 8: PLV')
% xlabel('Method')
ylabel('PLV','FontSize',13)
ylim([0,plv_max+0.02])
xlim([1,size(plv_MVMD,2)+size(plv_MSWD,2)])
start = (1 + length(freq_ranges))/2;
stop = 2*length(freq_ranges) + 1 - start;
ticks = start:length(freq_ranges):stop;
xticks(ticks)
xticklabels({'MVMD-PS','MSwD-PS'})
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',13)
xline((length(freq_ranges):length(freq_ranges):length(freq_ranges))+0.5,'--')
hold off
if save_plots
    exportgraphics(gcf,fullfile(pwd,'plots',['sim',num2str(casee),'_plv.png']),'Resolution',1000)
end