% Parameter sensitivity analysis for simulation 7

clear
clc
close all


%% Sim 7
save_plots = false;
casee = 7;
TR = 2;                                                                     % Repetition time
fs = 1/TR;                                                                  % Sampling frequency
f = 0.05;
t = 0:1/fs:668-1/fs;
N = 100;                                                 % number of repetition (realizations) 
SNR = -3;
min_freq = 0;

phi1 = 0;
phi2 = -pi/4;
phi3 = -pi/2;
phi4 = pi;
phi5 = 3*pi/4;
phi6 = pi/2;

delphis = [phi1-phi2,phi1-phi3,phi1-phi4,phi1-phi5,phi1-phi6,...
    phi2-phi3,phi2-phi4,phi2-phi5,phi2-phi6,phi3-phi4,...
    phi3-phi5,phi3-phi6,phi4-phi5,phi4-phi6,phi5-phi6];
COSDELPHI_grd = cos(delphis);

load(fullfile("simulated",['sim',num2str(casee),'.mat']))

indx = nchoosek(1:6,2);

%% MEMD
stop_vecs = {[0.075,0.75,0.075],[0.15,0.5,0.15],[0.2,0.9,0.2],[0.4,0.8,0.4],[0.3,0.3,0.3],[0.5,0.5,0.5]};
ndirs = [1,2,3,4,5,6,7,8];
if exist(fullfile("F:\MSWD_paper_new\parameter_sensitivity_analysis",['sim',num2str(casee),'_MEMD.mat']))
    load(fullfile("F:\MSWD_paper_new\parameter_sensitivity_analysis",['sim',num2str(casee),'_MEMD.mat']))
else
    imfs_MEMD = cell(length(ndirs),length(stop_vecs),N);
end
for k = 1:length(ndirs)
    for a = 1:length(stop_vecs)
        stop_vec = stop_vecs{a};
        ndir = ndirs(k);
        for m = 1:N
            if ~exist(fullfile("F:\MSWD_paper_new\parameter_sensitivity_analysis",['sim',num2str(casee),'_MEMD.mat']))
                stp_crit = 'stop';
                mode = 'na_snr';
                intensity_noise = 0.75; 
                n_channel_na = size(Data{m},2);  
                ndirr = ndir*n_channel_na;
                imfs = namemd(Data{m}, ndirr, stp_crit, stop_vec, mode, intensity_noise, n_channel_na);
                for i = 1:length(imfs)
                    imf(:,:,i) = imfs{i};
                end
                imfs_MEMD{k,a,m} = imf;
            else
                imf = imfs_MEMD{k,a,m};
            end
            [COSDELPHI,plv,mfreq] = phase_sync_analysis47(imf,"MVMD",indx,...
                fs,f,min_freq,1);
            clear imf
            for j = 1:size(indx,1)
                if ~all(isnan(COSDELPHI(:,j)))
                    rmse = sqrt(mean((COSDELPHI_grd(:,j) - COSDELPHI(:,j)).^2));
                    NRMSE_MEMD(k,a,m,j) = rmse;
                else
                    NRMSE_MEMD(k,a,m,j) = 1;
                end
            end
        end
        disp(['ndir:', num2str(ndirs(k)), ' stop_vec:', num2str(stop_vecs{a})])
    end
end

if ~exist(fullfile("F:\MSWD_paper_new\parameter_sensitivity_analysis",['sim',num2str(casee),'_MEMD.mat']))
    save(fullfile("F:\MSWD_paper_new\parameter_sensitivity_analysis",['sim',num2str(casee),'_MEMD.mat']),"imfs_MEMD")
end

NRMSE_MEMD = mean(mean(NRMSE_MEMD,4),3);


%% EWT
P_ths = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8];
freq_ress = [0.04,0.08,0.13,0.17,0.20,0.24];
if exist(fullfile("F:\MSWD_paper_new\parameter_sensitivity_analysis",['sim',num2str(casee),'_EWT.mat']))
    load(fullfile("F:\MSWD_paper_new\parameter_sensitivity_analysis",['sim',num2str(casee),'_EWT.mat']))
else
    imfs_EWT = cell(length(P_ths),length(freq_ress),N);
end
for k = 1:length(P_ths)
    for a = 1:length(freq_ress)
        P_th = P_ths(k);
        freq_res = freq_ress(a);
        for m = 1:N
            if ~exist(fullfile("F:\MSWD_paper_new\parameter_sensitivity_analysis",['sim',num2str(casee),'_EWT.mat']))
                imf = zeros(20,size(Data{m},1),size(Data{m},2));
                for i = 1:size(Data{m},2)
                    temp_imf = ewt(Data{m}(:,i),"PeakThresholdPercent",P_th,"FrequencyResolution",freq_res)';
                    mfreqs = zeros(1,size(temp_imf,1));
                    for j = 1:size(temp_imf,1)
                        mfreqs(j) = meanfreq(temp_imf(j,:),fs);
                    end
                    [~,pos] = min(abs(mfreqs - f));
                    temp = temp_imf(pos,:);
                    temp_imf(pos,:) = temp_imf(1,:);
                    temp_imf(1,:) = temp;
                    imf(1:size(temp_imf,1),:,i) = temp_imf;
                    clear temp_imf
                end
                imf(all(all(imf == 0,3),2),:,:) = [];
                imfs_EWT{k,a,m} = imf;
            else
                imf = imfs_EWT{k,a,m};
            end
            [COSDELPHI,plv,mfreq] = phase_sync_analysis47(imf,"MVMD",indx,...
                fs,f,min_freq,1);
            clear imf
            for j = 1:size(indx,1)
                if ~all(isnan(COSDELPHI(:,j)))
                    rmse = sqrt(mean((COSDELPHI_grd(:,j) - COSDELPHI(:,j)).^2));
                    NRMSE_EWT(k,a,m,j) = rmse;
                else
                    NRMSE_EWT(k,a,m,j) = 1;
                end
            end
        end
        disp(['P_th:', num2str(P_ths(k)), ' freq_res:', num2str(freq_ress(a))])
    end
end

if ~exist(fullfile("F:\MSWD_paper_new\parameter_sensitivity_analysis",['sim',num2str(casee),'_EWT.mat']))
    save(fullfile("F:\MSWD_paper_new\parameter_sensitivity_analysis",['sim',num2str(casee),'_EWT.mat']),"imfs_EWT")
end

NRMSE_EWT = mean(mean(NRMSE_EWT,4),3);

%% MVMD
Ks = 1:8;
alphas = [100,500,1000,2000,5000,10000];
if exist(fullfile("F:\MSWD_paper_new\parameter_sensitivity_analysis",['sim',num2str(casee),'_MVMD.mat']))
    load(fullfile("F:\MSWD_paper_new\parameter_sensitivity_analysis",['sim',num2str(casee),'_MVMD.mat']))
else
    imfs_MVMD = cell(length(Ks),length(alphas));
end
for k = 1:length(Ks)
    for a = 1:length(alphas)
        K = Ks(k);
        alpha = alphas(a);
        for m = 1:N
            if ~exist(fullfile("F:\MSWD_paper_new\parameter_sensitivity_analysis",['sim',num2str(casee),'_MVMD.mat']))
                tau = 0.01; DC = 1; init = 0; tol = 1e-9;
                imf = MVMD_new(Data{m},alpha,tau,K,DC,init,tol);
                imfs_MVMD{k,a,m} = imf;
            else
                imf = imfs_MVMD{k,a,m};
            end
            [COSDELPHI,plv,mfreq] = phase_sync_analysis47(imf,"MVMD",indx,...
                fs,f,min_freq,1);
            for j = 1:size(indx,1)
                if ~all(isnan(COSDELPHI(:,j)))
                    rmse = sqrt(mean((COSDELPHI_grd(:,j) - COSDELPHI(:,j)).^2));
                    NRMSE_MVMD(k,a,m,j) = rmse;
                else
                    NRMSE_MVMD(k,a,m,j) = 1;
                end
            end
        end
        disp(['K:', num2str(Ks(k)), ' alpha:', num2str(alphas(a))])
    end
end

if ~exist(fullfile("F:\MSWD_paper_new\parameter_sensitivity_analysis",['sim',num2str(casee),'_MVMD.mat']))
    save(fullfile("F:\MSWD_paper_new\parameter_sensitivity_analysis",['sim',num2str(casee),'_MVMD.mat']),"imfs_MVMD")
end

NRMSE_MVMD = mean(mean(NRMSE_MVMD,4),3);

%% MSWD
compStds = [0.002,0.005,0.01,0.05,0.1,0.15,0.2];
P_corr_imps = [0.01,0.02,0.05,0.07,0.1];

if exist(fullfile("F:\MSWD_paper_new\parameter_sensitivity_analysis",['sim',num2str(casee),'_MSWD.mat']))
    load(fullfile("F:\MSWD_paper_new\parameter_sensitivity_analysis",['sim',num2str(casee),'_MSWD.mat']))
else
    imfs_MSWD = cell(length(P_corr_imps),length(compStds));
end
for k = 1:length(P_corr_imps)
    for a = 1:length(compStds)
        wind = [];
        compStd = compStds(a);
        for m = 1:N
            if ~exist(fullfile("F:\MSWD_paper_new\parameter_sensitivity_analysis",['sim',num2str(casee),'_MSWD.mat']))
                p_value = 1e-5;
                param_struct  = struct('P_corr', 1, ...
                    'P_corr_imp',       P_corr_imps(k), ...
                    'StD_th',       compStd, ...
                    'Welch_window', wind, ...
                    'p_value',      1e-5);
                imf = MSWD(Data{m}, param_struct);
                imfs_MSWD{k,a,m} = imf;
            else
                imf = imfs_MSWD{k,a,m};
            end
            [COSDELPHI,plv,mfreq] = phase_sync_analysis47(imf,"MSWD",indx,...
                fs,f,min_freq,1);
            for j = 1:size(indx,1)
                if ~all(isnan(COSDELPHI(:,j)))
                    rmse = sqrt(mean((COSDELPHI_grd(:,j) - COSDELPHI(:,j)).^2));
                    NRMSE_MSWD(k,a,m,j) = rmse;
                else
                    NRMSE_MSWD(k,a,m,j) = 1;
                end
            end
        end
        disp(['P_corr_imp:', num2str(P_corr_imps(k)), ' compStd:', num2str(compStds(a))])
    end
end

if ~exist(fullfile("F:\MSWD_paper_new\parameter_sensitivity_analysis",['sim',num2str(casee),'_MSWD.mat']))
    save(fullfile("F:\MSWD_paper_new\parameter_sensitivity_analysis",['sim',num2str(casee),'_MSWD.mat']),"imfs_MSWD")
end

NRMSE_MSWD = mean(mean(NRMSE_MSWD,4),3);

%% MVMD-MSWD
maximum = max([max(NRMSE_MSWD,[],"all"),max(NRMSE_MVMD,[],"all")]);
minimum = min([min(NRMSE_MSWD,[],"all"),min(NRMSE_MVMD,[],"all")]);

figure
surf(Ks,1:length(alphas),NRMSE_MVMD');
ax = gca;
ax.FontSize = 12;
xlabel('K','FontSize',14)
xticks(Ks)
ylabel('alpha','FontSize',14)
yticks(1:length(alphas))
yticklabels({'100','500','1000','2000','5000','10000'})
zlim([minimum-0.2,maximum+0.2])
c=colorbar('Position',[1 0.001 0.001 0.01]); caxis([minimum maximum])
if save_plots
    exportgraphics(gcf,fullfile('plots_PSA_new',['MVMD_sim',num2str(casee),'_NRMSE.png']),'Resolution',1000)
end


figure
surf(1:length(P_corr_imps),1:length(compStds),NRMSE_MSWD');
ax = gca;
ax.FontSize = 12;
xlabel('Corr_t_h','FontSize',14)
xticks(1:length(P_corr_imps))
ylabel('StD_t_h','FontSize',14)
yticks(1:length(compStds))
yticklabels({'0.002','0.005','0.01','0.05','0.1','0.15','0.2'})
xticklabels({'0.01','0.02','0.05','0.07','0.1'})
zlim([minimum-0.2,maximum+0.2])
c=colorbar('Position',[1 0.001 0.001 0.01]); caxis([minimum maximum])
if save_plots
    exportgraphics(gcf,fullfile('plots_PSA_new',['MSWD_sim',num2str(casee),'_NRMSE.png']),'Resolution',1000)
end


%% MEMD-EWT-MVMD-MSWD
maximum = max([max(NRMSE_MEMD,[],"all"),max(NRMSE_EWT,[],"all"),max(NRMSE_MSWD,[],"all"),max(NRMSE_MVMD,[],"all")]);
minimum = min([min(NRMSE_MEMD,[],"all"),min(NRMSE_EWT,[],"all"),min(NRMSE_MSWD,[],"all"),min(NRMSE_MVMD,[],"all")]);

ndirs = size(Data{1},2)*ndirs;
figure
surf(ndirs,1:length(stop_vecs),NRMSE_MEMD');
ax = gca;
ax.FontSize = 12;
xlabel('$p$','Interpreter','latex','FontSize',18)
xticks(ndirs)
ylabel('[$sd_1$,$sd_2$,$tol$]','Interpreter','latex','FontSize',18)
yticks(1:length(stop_vecs))
zlim([minimum-0.2,maximum+0.2])
yticklabels({'[0.075,0.75,0.075]','[0.15,0.5,0.15]','[0.2,0.9,0.2]','[0.4,0.8,0.4]','[0.3,0.3,0.3]','[0.5,0.5,0.5]'})
c=colorbar('Position',[1 0.001 0.001 0.01]); caxis([minimum maximum])
if save_plots
    exportgraphics(gcf,fullfile(pwd,'plots_PSA_new',['MEMD_sim',num2str(casee),'_NRMSE.png']),'Resolution',1000)
end

figure
surf(1:length(P_ths),1:length(freq_ress),NRMSE_EWT');
ax = gca;
ax.FontSize = 12;
xlabel('P_t_h','FontSize',14)
xticks(1:length(P_ths))
xticklabels({'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8'})
ylabel('freqRes','FontSize',14)
yticks(1:length(freq_ress))
zlim([minimum-0.2,maximum+0.2])
yticklabels({'0.04','0.08','0.13','0.17','0.20','0.25'})
c=colorbar('Position',[1 0.001 0.001 0.01]); caxis([minimum maximum])
if save_plots
    exportgraphics(gcf,fullfile(pwd,'plots_PSA_new',['EWT_sim',num2str(casee),'_NRMSE.png']),'Resolution',1000)
end

figure
surf(Ks,1:length(alphas),NRMSE_MVMD');
ax = gca;
ax.FontSize = 12;
xlabel('K','FontSize',14)
xticks(Ks)
ylabel('alpha','FontSize',14)
yticks(1:length(alphas))
zlim([minimum-0.2,maximum+0.2])
yticklabels({'100','500','1000','2000','5000','10000'})
c=colorbar('Position',[1 0.001 0.001 0.01]); caxis([minimum maximum])
if save_plots
    exportgraphics(gcf,fullfile(pwd,'plots_PSA_new',['MVMD_all_sim',num2str(casee),'_NRMSE.png']),'Resolution',1000)
end

figure
surf(1:length(P_corr_imps),1:length(compStds),NRMSE_MSWD');
ax = gca;
ax.FontSize = 12; 
xlabel('Corr_t_h','FontSize',14)
xticks(1:length(P_corr_imps))
ylabel('StD_t_h','FontSize',14)
yticks(1:length(compStds))
zlim([minimum-0.2,maximum+0.2])
yticklabels({'0.002','0.005','0.01','0.05','0.1','0.15','0.2'})
xticklabels({'0.01','0.02','0.05','0.07','0.1'})
c=colorbar('Position',[1 0.001 0.001 0.01]); caxis([minimum maximum])
if save_plots
    exportgraphics(gcf,fullfile(pwd,'plots_PSA_new',['MSWD_all_sim',num2str(casee),'_NRMSE.png']),'Resolution',1000)
end

