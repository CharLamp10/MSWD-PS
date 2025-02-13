clear
clc
close all


%% Sim 3
save_plots = 1;
casee = 3;
TR = 2;                                                   % Repetition time 
fs = 1/TR;                                                % Sampling frequency
t = 0:1/fs:668-1/fs;
f = 0.05;
N = 100;                                                 % number of repetition (realizations)
w = [30 60 120];  
freq_ranges = [0.01,0.03,0.05,0.055,0.07,0.09];
min_freq = 0;                                        % window sizes for the Windowed Phase Sync. Measures

true_plv = zeros(1,length(freq_ranges));
delphi = 2*pi./(1+exp(-0.01*(t-334)));

x = cos(2*pi*f*t);                                          % first signal  (pure tone - ie monocomponent)
y = cos(2*pi*f*t + delphi) + cos(2*pi*f*1.1*t + delphi);    % second signal (multicomponent)
y_plv = cos(2*pi*f*t + delphi);
pwelch_window = length(x);
phase_x = angle(hilbert(x)); phase_y = angle(hilbert(y_plv));
dPhi = phase_x - phase_y;
e = exp(1i*dPhi);
pos = find(freq_ranges == f);
true_plv(pos) = abs(mean(e));
true_cosdelphi = cos(phase_x - phase_y);
filename = fullfile(pwd,'plots_new',['sim',num2str(casee) , '_ground_truth.png']);
simulation2_groundTruthPlot(true_cosdelphi,t)

load(fullfile("decomposed",['sim',num2str(casee),'_MEMD.mat']))
load(fullfile("decomposed",['sim',num2str(casee),'_EWT.mat']))
load(fullfile("decomposed",['sim',num2str(casee),'_RRCNN.mat']))
load(fullfile("decomposed",['sim',num2str(casee),'_MVMD.mat']))
load(fullfile("decomposed",['sim',num2str(casee),'_MSWD.mat']))

%% MEMD
for m = 1:N
    imf = imfs_MEMD{m};
    if m == 1
        plv_length = length(freq_ranges);
    else
        plv_length = length(plv(m-1,:));
    end
    [COSDELPHI11,plv1,mfreq1] = phase_sync_analysis12356(imf,...
    "MVMD",fs,min_freq,plv_length,f);
    COSDELPHI1(:,m) = COSDELPHI11;
    mfreq{m} = mfreq1;
    plv(m,1:length(plv1)) = plv1;
end
[plv_MEMD,counts_MEMD] = sort_values(mfreq,freq_ranges,plv);
[mean_plv_MEMD,lower_plv_MEMD,upper_plv_MEMD] = confidence_interval(plv_MEMD);

for m = 1:N
    cosdelphi_memd(m,:) = COSDELPHI1(:,m);
end
disp('MEMD')
clear COSDELPHI1 CCORSW mfreq plv

%% EWT
for m = 1:N
    imf = imfs_EWT{m};
    if m == 1
        plv_length = length(freq_ranges);
    else
        plv_length = length(plv(m-1,:));
    end
    [COSDELPHI11,plv1,mfreq1] = phase_sync_analysis12356(imf,...
    "MVMD",fs,min_freq,plv_length,f);
    COSDELPHI1(:,m) = COSDELPHI11;
    mfreq{m} = mfreq1;
    plv(m,1:length(plv1)) = plv1;
end
[plv_EWT,counts_EWT] = sort_values(mfreq,freq_ranges,plv);
[mean_plv_EWT,lower_plv_EWT,upper_plv_EWT] = confidence_interval(plv_EWT);

for m = 1:N
    cosdelphi_ewt(m,:) = COSDELPHI1(:,m);
end
disp('EWT')
clear COSDELPHI1 CCORSW mfreq plv


%% RRCNN
for m = 1:N
    imf = imfs_RRCNN{m};
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
    if m == 1
        plv_length = length(freq_ranges);
    else
        plv_length = length(plv(m-1,:));
    end
    [COSDELPHI11,plv1,mfreq1] = phase_sync_analysis12356(imf,...
    "MVMD",fs,min_freq,plv_length,f);
    COSDELPHI1(:,m) = COSDELPHI11;
    mfreq{m} = mfreq1;
    plv(m,1:length(plv1)) = plv1;
end
[plv_RRCNN,counts_RRCNN] = sort_values(mfreq,freq_ranges,plv);
[mean_plv_RRCNN,lower_plv_RRCNN,upper_plv_RRCNN] = confidence_interval(plv_RRCNN);

for m = 1:N
    cosdelphi_rrcnn(m,:) = COSDELPHI1(:,m);
end
disp('RRCNN')
clear COSDELPHI1 CCORSW mfreq plv

%% MVMD
for m = 1:N
    imf = imfs_MVMD{m};
    if m == 1
        plv_length = length(freq_ranges);
    else
        plv_length = length(plv(m-1,:));
    end
    [COSDELPHI11,plv1,mfreq1] = phase_sync_analysis12356(imf,...
    "MVMD",fs,min_freq,plv_length,f);
    COSDELPHI1(:,m) = COSDELPHI11;
    mfreq{m} = mfreq1;
    plv(m,1:length(plv1)) = plv1;
end
[plv_MVMD,counts_MVMD] = sort_values(mfreq,freq_ranges,plv);
[mean_plv_MVMD,lower_plv_MVMD,upper_plv_MVMD] = confidence_interval(plv_MVMD);

for m = 1:N
    cosdelphi_mvmd(m,:) = COSDELPHI1(:,m);
end
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
    "MSWD",fs,min_freq,plv_length,f);
    COSDELPHI1(:,m) = COSDELPHI11;
    mfreq{m} = mfreq1;
    plv(m,1:length(plv1)) = plv1;
end
[plv_MSWD,counts_MSWD] = sort_values(mfreq,freq_ranges,plv);
[mean_plv_MSWD,lower_plv_MSWD,upper_plv_MSWD] = confidence_interval(plv_MSWD);

for m = 1:N
    cosdelphi_mswd(m,:) = COSDELPHI1(:,m);
end
disp('MSWD')

%% PLOTS
filename = fullfile(pwd,'plots_new',['sim', num2str(casee), '_TVPS.png']);
if save_plots
    simulation12356Plots(cosdelphi_memd,cosdelphi_ewt,cosdelphi_rrcnn,cosdelphi_mvmd,cosdelphi_mswd,t,filename)
else
    simulation12356Plots(cosdelphi_memd,cosdelphi_ewt,cosdelphi_rrcnn,cosdelphi_mvmd,cosdelphi_mswd,t)
end



E_plv_MEMD = mean_plv_MEMD - lower_plv_MEMD;
E_plv_EWT = mean_plv_EWT - lower_plv_EWT;
E_plv_RRCNN = mean_plv_RRCNN - lower_plv_RRCNN;
E_plv_MVMD = mean_plv_MVMD - lower_plv_MVMD;
E_plv_MSWD = mean_plv_MSWD - lower_plv_MSWD;


plv_min = min([min(lower_plv_MEMD),min(lower_plv_EWT),min(lower_plv_RRCNN),min(lower_plv_MVMD),min(lower_plv_MSWD),min(true_plv)]);
plv_max = max([max(lower_plv_MEMD),max(lower_plv_EWT),max(lower_plv_RRCNN),max(upper_plv_MVMD),max(upper_plv_MSWD),max(true_plv)]);


figure
errorbar(1:size(plv_MEMD,2),mean_plv_MEMD,E_plv_MEMD,'mo',"LineStyle","none",'LineWidth',1)
hold on
errorbar(size(plv_MEMD,2)+1:2*size(plv_MEMD,2),mean_plv_EWT,E_plv_EWT,'co',"LineStyle","none",'LineWidth',1)
errorbar(2*size(plv_MEMD,2)+1:3*size(plv_MEMD,2),mean_plv_RRCNN,E_plv_RRCNN,'yo',"LineStyle","none",'LineWidth',1)
errorbar(3*size(plv_MEMD,2)+1:4*size(plv_MEMD,2),mean_plv_MVMD,E_plv_MVMD,'bo',"LineStyle","none",'LineWidth',1)
errorbar(4*size(plv_MEMD,2)+1:5*size(plv_MEMD,2),mean_plv_MSWD,E_plv_MSWD,'go',"LineStyle","none",'LineWidth',1)
plot(1:5*size(plv_MEMD,2),[true_plv,true_plv,true_plv,true_plv,true_plv],'ro','LineWidth',1)
ylabel('PLV','FontSize',16)
ylim([0,plv_max+0.02])
xlim([1,5*size(plv_MEMD,2)])
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





