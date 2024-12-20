% Parameter sensitivity analysis for simulation 4

clear
clc
close all


%% Sim 4
save_plots = 0;
casee = 4;
TR = 2;                                                                     % Repetition time
fs = 1/TR;                                                                  % Sampling frequency
f = 0.05;
t = 0:1/fs:1000-1/fs;
N = 100;                                                 % number of repetition (realizations) 
nS = 3;
SNR = -3;
min_freq = 0;

phi1 = pi*((t-100>0).*(t-250<0) + (t-300>0).*(t-500<0) + (t-600>0).*(t-800<0));
phi2 = pi*((t-100>0).*(t-250<0) +                      - (t-600>0).*(t-800<0));
phi3 = pi*(                     - (t-300>0).*(t-500<0) - (t-600>0).*(t-800<0));
COSDELPHI_grd(:,1) = cos(phi1 - phi2);
COSDELPHI_grd(:,2) = cos(phi1 - phi3);
COSDELPHI_grd(:,3) = cos(phi2 - phi3);

load(fullfile("simulated",['sim',num2str(casee),'.mat']))

indx = nchoosek(1:size(Data{1},2),2);

%% MVMD
Ks = 1:8;
alphas = [100,500,1000,2000,5000,10000];
if exist(fullfile("E:\MSWD_paper_new\parameter_sensitivity_analysis",['sim',num2str(casee),'_MVMD.mat']))
    load(fullfile("E:\MSWD_paper_new\parameter_sensitivity_analysis",['sim',num2str(casee),'_MVMD.mat']))
else
    imfs_MVMD = cell(length(Ks),length(alphas));
end
for k = 1:length(Ks)
    for a = 1:length(alphas)
        K = Ks(k);
        alpha = alphas(a);
        for m = 1:N
            if ~exist(fullfile("E:\MSWD_paper_new\parameter_sensitivity_analysis",['sim',num2str(casee),'_MVMD.mat']))
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

if ~exist(fullfile("E:\MSWD_paper_new\parameter_sensitivity_analysis",['sim',num2str(casee),'_MVMD.mat']))
    save(fullfile("E:\MSWD_paper_new\parameter_sensitivity_analysis",['sim',num2str(casee),'_MVMD.mat']),"imfs_MVMD")
end

NRMSE_MVMD = mean(mean(NRMSE_MVMD,4),3);


%% MSWD
compStds = [0.002,0.005,0.01,0.05,0.1,0.15,0.2];
P_corr_imps = [0.01,0.02,0.05,0.07,0.1];
if exist(fullfile("E:\MSWD_paper_new\parameter_sensitivity_analysis",['sim',num2str(casee),'_MSWD.mat']))
    load(fullfile("E:\MSWD_paper_new\parameter_sensitivity_analysis",['sim',num2str(casee),'_MSWD.mat']))
else
    imfs_MSWD = cell(length(P_corr_imps),length(compStds));
end
for k = 1:length(P_corr_imps)
    for a = 1:length(compStds)
        wind = [];
        compStd = compStds(a);
        for m = 1:N
            if ~exist(fullfile("E:\MSWD_paper_new\parameter_sensitivity_analysis",['sim',num2str(casee),'_MSWD.mat']))
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

if ~exist(fullfile("E:\MSWD_paper_new\parameter_sensitivity_analysis",['sim',num2str(casee),'_MSWD.mat']))
    save(fullfile("E:\MSWD_paper_new\parameter_sensitivity_analysis",['sim',num2str(casee),'_MSWD.mat']),"imfs_MSWD")
end

NRMSE_MSWD = mean(mean(NRMSE_MSWD,4),3);

maximum = max([max(NRMSE_MSWD,[],"all"),max(NRMSE_MVMD,[],"all")]) + 0.2;
minimum = min([min(NRMSE_MSWD,[],"all"),min(NRMSE_MVMD,[],"all")]) - 0.2;

figure
surf(Ks,1:length(alphas),NRMSE_MVMD');
ax = gca;
ax.FontSize = 12;
xlabel('K','FontSize',14)
xticks(Ks)
ylabel('alpha','FontSize',14)
yticks(1:length(alphas))
yticklabels({'100','500','1000','2000','5000','10000','20000'})
zlim([minimum,maximum])
c=colorbar('Position',[1 0.001 0.001 0.01]); caxis([minimum+0.2 maximum-0.2])
if save_plots
    exportgraphics(gcf,fullfile(pwd,'plots_PSA',['MVMD_sim',num2str(casee),'_NRMSE.png']),'Resolution',1000)
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
zlim([minimum,maximum])
c=colorbar('Position',[1 0.001 0.001 0.01]); caxis([minimum+0.2 maximum-0.2])
if save_plots
    exportgraphics(gcf,fullfile(pwd,'plots_PSA',['MSWD_sim',num2str(casee),'_NRMSE.png']),'Resolution',1000)
end



