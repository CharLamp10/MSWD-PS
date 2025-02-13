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
% filename = fullfile(pwd,'plots_new',['sim',num2str(casee) , '_ground_truth.png']);
% simulation2_groundTruthPlot(true_cosdelphi,t,filename)

load(fullfile("decomposed",['sim',num2str(casee),'_MEMD.mat']))
load(fullfile("decomposed",['sim',num2str(casee),'_EWT.mat']))
load(fullfile("decomposed",['sim',num2str(casee),'_IRCNN.mat']))
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
    "MVMD",fs,min_freq,plv_length,f1);
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
    "MVMD",fs,min_freq,plv_length,f1);
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


%% IRCNN
for m = 1:N
    imf = imfs_IRCNN{m};
    for i = 1:size(imf,3)
        temp_imf = imf(:,:,i);
        mfreqs = zeros(1,size(temp_imf,1));
        for j = 1:size(temp_imf,1)
            mfreqs(j) = meanfreq(temp_imf(j,:),fs);
        end
        [val1,pos1] = min(abs(mfreqs - f1));
        [val2,pos2] = min(abs(mfreqs - f2));
        if pos1 ~= pos2
            temp1 = temp_imf(pos1,:);
            temp2 = temp_imf(pos2,:);
        else
            if val1 < val2
                temp1 = temp_imf(pos1,:);
                mfreqs(pos1) = 1e10;
                [~,pos2] = min(abs(mfreqs - f2));
                temp2 = temp_imf(pos2,:);
            else
                temp2 = temp_imf(pos2,:);
                mfreqs(pos2) = 1e10;
                [~,pos1] = min(abs(mfreqs - f2));
                temp1 = temp_imf(pos1,:);
            end
        end
        temp_imf(pos1,:) = temp_imf(1,:);
        temp_imf(pos2,:) = temp_imf(2,:);
        temp_imf(1,:) = temp1;
        temp_imf(2,:) = temp2;
        imf(1:size(temp_imf,1),1:size(temp_imf,2),i) = temp_imf;
        clear temp_imf
    end
    if m == 1
        plv_length = length(freq_ranges);
    else
        plv_length = length(plv(m-1,:));
    end
    [COSDELPHI11,plv1,mfreq1] = phase_sync_analysis12356(imf,...
    "MVMD",fs,min_freq,plv_length,f1);
    COSDELPHI1(:,m) = COSDELPHI11;
    mfreq{m} = mfreq1;
    plv(m,1:length(plv1)) = plv1;
end
[plv_IRCNN,counts_IRCNN] = sort_values(mfreq,freq_ranges,plv);
[mean_plv_IRCNN,lower_plv_IRCNN,upper_plv_IRCNN] = confidence_interval(plv_IRCNN);

for m = 1:N
    cosdelphi_ircnn(m,:) = COSDELPHI1(:,m);
end
disp('IRCNN')
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
filename = fullfile(pwd,'plots_new',['sim', num2str(casee), '_TVPS.png']);
if save_plots
    simulation12356Plots(cosdelphi_memd,cosdelphi_ewt,cosdelphi_ircnn,cosdelphi_mvmd,cosdelphi_mswd,t,filename)
else
    simulation12356Plots(cosdelphi_memd,cosdelphi_ewt,cosdelphi_ircnn,cosdelphi_mvmd,cosdelphi_mswd,t)
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
errorbar(1:size(plv_MEMD,2),mean_plv_MEMD,E_plv_MEMD,'mo',"LineStyle","none",'LineWidth',1)
hold on
errorbar(size(plv_MEMD,2)+1:2*size(plv_MEMD,2),mean_plv_EWT,E_plv_EWT,'co',"LineStyle","none",'LineWidth',1)
errorbar(2*size(plv_MEMD,2)+1:3*size(plv_MEMD,2),mean_plv_IRCNN,E_plv_IRCNN,'yo',"LineStyle","none",'LineWidth',1)
errorbar(3*size(plv_MEMD,2)+1:4*size(plv_MEMD,2),mean_plv_MVMD,E_plv_MVMD,'bo',"LineStyle","none",'LineWidth',1)
errorbar(4*size(plv_MEMD,2)+1:5*size(plv_MEMD,2),mean_plv_MSWD,E_plv_MSWD,'go',"LineStyle","none",'LineWidth',1)
plot([true_plv,true_plv,true_plv,true_plv,true_plv],'ro','LineWidth',1)
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