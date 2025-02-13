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

load(fullfile("decomposed",['sim',num2str(casee),'_MEMD.mat']))
load(fullfile("decomposed",['sim',num2str(casee),'_EWT.mat']))
load(fullfile("decomposed",['sim',num2str(casee),'_IRCNN.mat']))
load(fullfile("decomposed",['sim',num2str(casee),'_MVMD.mat']))
load(fullfile("decomposed",['sim',num2str(casee),'_MSWD.mat']))


%% MEMD
for m = 1:N
    imf = imfs_MEMD{m};
    imf = squeeze(imf);
    if m == 1
        plv_length = [length(freq_ranges),size(indx,1)];
    else
        plv_length = [length(plv(m-1,:,1)),size(indx,1)];
    end
    [COSDELPHI,plv1,mfreq1] = phase_sync_analysis47(imf,"MVMD",indx,...
     fs,f,min_freq,plv_length);
    COSDELPHI_MEMD{m} = COSDELPHI;
    for i = 1:size(COSDELPHI,2)
        nrmse_MEMD(m,i) = rmse(COSDELPHI(:,i),cosdelphi_grd(:,i));
    end
    plv(m,1:size(plv1,1),1:size(plv1,2)) = plv1;
    mfreq{m} = mfreq1;
end

for i = 1:size(indx,1)
    [plv_MEMD(:,:,i),counts_MEMD(:,:,i)] = ...
        sort_values(mfreq,freq_ranges,plv(:,:,i));
end  
plv_MEMD = mean(plv_MEMD,3);
counts_MEMD = mean(counts_MEMD,3);
[mean_plv_MEMD,lower_plv_MEMD,upper_plv_MEMD] = confidence_interval(plv_MEMD);

filename1 = fullfile(pwd,'plots_new',['sim', num2str(casee), '_MEMD']);
if save_plots
    simulation7Plots(COSDELPHI_MEMD,cosdelphi_grd,[],[],[],[],[],t,filename1);
else
    simulation7Plots(COSDELPHI_MEMD,cosdelphi_grd,[],[],[],[],[],t);
end
disp('MEMD')
clear mfreq plv

%% EWT
for m = 1:N
    imf = imfs_EWT{m};
    imf = squeeze(imf);
    if m == 1
        plv_length = [length(freq_ranges),size(indx,1)];
    else
        plv_length = [length(plv(m-1,:,1)),size(indx,1)];
    end
    [COSDELPHI,plv1,mfreq1] = phase_sync_analysis47(imf,"MVMD",indx,...
     fs,f,min_freq,plv_length);
    COSDELPHI_EWT{m} = COSDELPHI;
    for i = 1:size(COSDELPHI,2)
        nrmse_EWT(m,i) = rmse(COSDELPHI(:,i),cosdelphi_grd(:,i));
    end
    plv(m,1:size(plv1,1),1:size(plv1,2)) = plv1;
    mfreq{m} = mfreq1;
end

for i = 1:size(indx,1)
    [plv_EWT(:,:,i),counts_EWT(:,:,i)] = ...
        sort_values(mfreq,freq_ranges,plv(:,:,i));
end  
plv_EWT = mean(plv_EWT,3);
counts_EWT = mean(counts_EWT,3);
[mean_plv_EWT,lower_plv_EWT,upper_plv_EWT] = confidence_interval(plv_EWT);

filename1 = fullfile(pwd,'plots_new',['sim', num2str(casee), '_EWT']);
if save_plots
    simulation7Plots(COSDELPHI_EWT,cosdelphi_grd,[],[],[],[],[],t,filename1);
else
    simulation7Plots(COSDELPHI_EWT,cosdelphi_grd,[],[],[],[],[],t);
end
disp('EWT')
clear mfreq plv

%% IRCNN
for m = 1:N
    imf = imfs_IRCNN{m};
    for i = 1:size(imf,3)
        temp_imf = imf(:,:,i);
        mfreqs = zeros(1,size(temp_imf,1));
        for j = 1:size(temp_imf,1)
            mfreqs(j) = meanfreq(temp_imf(j,:),fs);
        end
        [~,pos] = min(abs(mfreqs - f));
        temp = temp_imf(pos,:);
        temp_imf(pos,:) = temp_imf(1,:);
        temp_imf(1,:) = temp;
        imf(:,:,i) = temp_imf;
        clear temp_imf
    end
    imf = squeeze(imf);
    if m == 1
        plv_length = [length(freq_ranges),size(indx,1)];
    else
        plv_length = [length(plv(m-1,:,1)),size(indx,1)];
    end
    [COSDELPHI,plv1,mfreq1] = phase_sync_analysis47(imf,"MVMD",indx,...
     fs,f,min_freq,plv_length);
    COSDELPHI_IRCNN{m} = COSDELPHI;
    for i = 1:size(COSDELPHI,2)
        nrmse_IRCNN(m,i) = rmse(COSDELPHI(:,i),cosdelphi_grd(:,i));
    end
    plv(m,1:size(plv1,1),1:size(plv1,2)) = plv1;
    mfreq{m} = mfreq1;
end

for i = 1:size(indx,1)
    [plv_IRCNN(:,:,i),counts_IRCNN(:,:,i)] = ...
        sort_values(mfreq,freq_ranges,plv(:,:,i));
end  
plv_IRCNN = mean(plv_IRCNN,3);
counts_IRCNN = mean(counts_IRCNN,3);
[mean_plv_IRCNN,lower_plv_IRCNN,upper_plv_IRCNN] = confidence_interval(plv_IRCNN);

filename1 = fullfile(pwd,'plots_new',['sim', num2str(casee), '_IRCNN']);
if save_plots
    simulation7Plots(COSDELPHI_IRCNN,cosdelphi_grd,[],[],[],[],[],t,filename1);
else
    simulation7Plots(COSDELPHI_IRCNN,cosdelphi_grd,[],[],[],[],[],t);
end
disp('EWT')
clear mfreq plv

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

filename1 = fullfile(pwd,'plots_new',['sim', num2str(casee), '_MVMD']);
if save_plots
    simulation7Plots(COSDELPHI_MVMD,cosdelphi_grd,[],[],[],[],[],t,filename1);
else
    simulation7Plots(COSDELPHI_MVMD,cosdelphi_grd,[],[],[],[],[],t);
end
disp('MVMD')
clear mfreq plv


%% MSWD
for m = 1:N
    imf = imfs_MSWD{m};
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


%% PLOTS
filename1 = fullfile(pwd,'plots_new',['sim', num2str(casee), '_MSWD']);
filename2 = fullfile(pwd,'plots_new',['sim', num2str(casee), '_RMSE']);
if save_plots
    simulation7Plots(COSDELPHI_MSWD,cosdelphi_grd,mean(nrmse_MEMD),mean(nrmse_EWT),mean(nrmse_IRCNN),mean(nrmse_MVMD),mean(nrmse_MSWD),t,filename1,filename2);
else
    simulation7Plots(COSDELPHI_MSWD,cosdelphi_grd,mean(nrmse_MEMD),mean(nrmse_EWT),mean(nrmse_IRCNN),mean(nrmse_MVMD),mean(nrmse_MSWD),t);
end
disp('MSWD')


E_plv_MEMD = mean_plv_MEMD - lower_plv_MEMD;
E_plv_EWT = mean_plv_EWT - lower_plv_EWT;
E_plv_IRCNN = mean_plv_IRCNN - lower_plv_IRCNN;
E_plv_MVMD = mean_plv_MVMD - lower_plv_MVMD;
E_plv_MSWD = mean_plv_MSWD - lower_plv_MSWD;

plv_min = min([min(lower_plv_MEMD),min(lower_plv_EWT),min(lower_plv_IRCNN) ...
    min(lower_plv_MVMD),min(lower_plv_MSWD),min(true_plv)]);
plv_max = max([max(lower_plv_MEMD),max(lower_plv_EWT),max(lower_plv_IRCNN) ...
    max(upper_plv_MVMD),max(upper_plv_MSWD),max(true_plv)]);

figure
errorbar(1:size(plv_MEMD,2),mean_plv_MEMD,E_plv_MEMD,'mo',"LineStyle","none")
hold on
errorbar(size(plv_MEMD,2)+1:2*size(plv_MEMD,2),mean_plv_EWT,E_plv_EWT,'co',"LineStyle","none")
errorbar(2*size(plv_MEMD,2)+1:3*size(plv_MEMD,2),mean_plv_IRCNN,E_plv_IRCNN,'yo',"LineStyle","none")
errorbar(3*size(plv_MEMD,2)+1:4*size(plv_MEMD,2),mean_plv_MVMD,E_plv_MVMD,'bo',"LineStyle","none")
errorbar(4*size(plv_MEMD,2)+1:5*size(plv_MEMD,2),mean_plv_MSWD,E_plv_MSWD,'go',"LineStyle","none")
plot([true_plv,true_plv,true_plv,true_plv,true_plv],'ro')
ylabel('PLV','FontSize',13)
ylim([0,plv_max+0.02])
xlim([1,5*size(plv_MVMD,2)])
start = (1 + length(freq_ranges))/2;
stop = 5*length(freq_ranges) + 1 - start;
ticks = start:length(freq_ranges):stop;
xticks(ticks)
xticklabels({'MEMD-PS','EWT-PS','IRCNN-PS','MVMD-PS','MSwD-PS'})
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',11,'FontWeight', 'bold')
xline((length(freq_ranges):length(freq_ranges):5*length(freq_ranges))+0.5,'--')
hold off
if save_plots
    exportgraphics(gcf,fullfile(pwd,'plots_new',['sim',num2str(casee),'_plv.png']),'Resolution',1000)
end