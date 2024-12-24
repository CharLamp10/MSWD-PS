% In this script, fMRI data from the HCP project are decomposed using MVMD
% and MSWD and CRP is calculated for dominant component

clear
clc

save_res = 1;
numICs = 100;
visit = 'rest1';
path_data = 'C:\Users\Hp\Desktop\MSWD_paper_files';
path_save = 'C:\Users\Hp\Desktop\MSWD_paper_files';

path_tmaps = fullfile(path_data,'dyn_tmap_comps_centered.nii');
path_hcp_average = fullfile(path_data,'HCP_PTN1200\groupICA\groupICA_3T_HCP1200_MSMAll_d100.ica\melodic_IC_sum.nii');
path_data = fullfile(path_data,'HCP_PTN1200\NodeTimeseries_3T_HCP1200_MSMAll_ICAd100_ts2\node_timeseries\3T_HCP1200_MSMAll_d100_ts2');


tmaps = niftiread(path_tmaps);
hcp_average = niftiread(path_hcp_average);

SC = [57,8,13,15];
AUD = [99,98];
SM = [2,10,14,58,37,45,30,9];
VIS = [77,79,60,46,27,43,82,89,80,61];
CC = [67,93,56,90,71,76,47,88,21,59,50,29,66,81];
DM = [20,44,39,28,64,26,75,48,87];
CB = [12,7,32];

[~,~,~,nICs] = size(hcp_average);
[~,~,~,nICs1] = size(tmaps);

cov_term = zeros(nICs,nICs1);
for i = 1:nICs
    resizedVolume = imresize3(hcp_average(:,:,:,i), size(tmaps(:,:,:,1)));
    var_hcp = var(resizedVolume(:));
    for j = 1:nICs1
        tmap = tmaps(:,:,:,j);
        covariance_matrix = corrcoef(resizedVolume(:), tmap(:));
        cov_term(i,j) = abs(covariance_matrix(1, 2));%^2/var_hcp;
    end
end

cov_term_new_max = zeros(nICs,7);
cov_term_new_max(:,1) = max(cov_term(:,SC),[],2);
cov_term_new_max(:,2) = max(cov_term(:,AUD),[],2);
cov_term_new_max(:,3) = max(cov_term(:,SM),[],2);
cov_term_new_max(:,4) = max(cov_term(:,VIS),[],2);
cov_term_new_max(:,5) = max(cov_term(:,CC),[],2);
cov_term_new_max(:,6) = max(cov_term(:,DM),[],2);
cov_term_new_max(:,7) = max(cov_term(:,CB),[],2);

cov_term_new_sum = zeros(nICs,7);
cov_term_new_sum(:,1) = sum(maxk(cov_term(:,SC)',2));
cov_term_new_sum(:,2) = sum(maxk(cov_term(:,AUD)',2));
cov_term_new_sum(:,3) = sum(maxk(cov_term(:,SM)',2));
cov_term_new_sum(:,4) = sum(maxk(cov_term(:,VIS)',2));
cov_term_new_sum(:,5) = sum(maxk(cov_term(:,CC)',2));
cov_term_new_sum(:,6) = sum(maxk(cov_term(:,DM)',2));
cov_term_new_sum(:,7) = sum(maxk(cov_term(:,CB)',2));

maximums = max(cov_term_new_max,[],2);
for i = 1:size(cov_term_new_max,1)
    positions(i) = find(cov_term_new_max(i,:) == maximums(i));
end
[positions_sort,I] = sort(positions);

dir_data = dir(path_data);

fs = 1/0.72;
thresh_same_imf = 10;
min_freq = 0.01;

numSubjects = 98;
compStds = zeros(1,numSubjects);
winds = zeros(1,numSubjects);
minPeaks = zeros(1,numSubjects);
imfs_MSWD = cell(1,numSubjects);
imfs_MVMD = cell(1,numSubjects);

if exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\MSWD_HCP_PSA.mat']))
    resMSWD = load(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\MSWD_HCP_PSA.mat']));
end
if exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\MVMD_HCP_PSA.mat']))
    resMVMD = load(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\MVMD_HCP_PSA.mat']));
end

load("names_98_subjects.mat");


inds_MSWD = [1]; %In MSWD we get always the first component

indx = nchoosek(1:numICs,2);

Ks = 1:2:7;
alphas = [500,2000,5000];
compStds = [0.002,0.01,0.1,0.2];
P_corr_imps = [0.01,0.07,0.1];

% for k = 1:length(compStds)
%     for a = 1:length(P_corr_imps)
%         for i = 1:numSubjects
%             name = names{i};
%             data = load(fullfile(path_data,name));
%             %% MSWD
%             if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\MSWD_HCP_PSA.mat']))
%                 data_MSWD = data(1:1200,I);
%                 data_MSWD = data_MSWD./max(abs(data_MSWD));
%                 p_value = 1e-5;
%                 P_corr = 1;
%                 P_corr_imp = P_corr_imps(a);
%                 wind = [];
%                 compStd = compStds(k);
%                 param_struct  = struct('P_corr', P_corr, ...
%                     'P_corr_imp',   P_corr_imp,...
%                     'StD_th',       compStd, ...
%                     'Welch_window', wind, ...
%                     'p_value',      1e-5);
%                 imf = MSWD(data_MSWD, param_struct);
%                 imfs_MSWD{i,k,a} = imf;
%                 disp(['P_corr_imp:', num2str(P_corr_imps(a)), ' compStd:', num2str(compStds(k))])
%             elseif ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_MSWD_PSA.mat']))
%                 imf = resMSWD.imfs_MSWD{i,k,a};
%             end
%             if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_MSWD_PSA.mat'])) && i <= 50
%                 COSDELPHI = phase_sync_analysis_HCP(imf,'MSWD',indx,...
%                      inds_MSWD);
%                 COSDELPHI1_MSWD{i,k,a} = COSDELPHI;
%             end
%         end
%     end
% end

for k = 1:length(Ks)
    for a = 1:length(alphas)
        for i = 1:numSubjects
            name = names{i};
            data = load(fullfile(path_data,name));
            %% MVMD
            if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\MVMD_HCP_PSA.mat']))
                data_MVMD = data(1:1200,I);
                data_MVMD = data_MVMD./max(abs(data_MVMD));
                tau = 0; DC = 1; init = 0; tol = 1e-9;
                K = Ks(k); alpha = alphas(a);
                imf = MVMD_new(data_MVMD,alpha,tau,K,DC,init,tol);
                imfs_MVMD{i,k,a} = imf;
                disp(['K:', num2str(Ks(k)), ' alpha:', num2str(alphas(a))])
            end
        end
    end
end
if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\MVMD_HCP_PSA.mat']))
    resMVMD.imfs_MVMD = imfs_MVMD;
end
if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_MVMD_PSA.mat']))
    for k = 1:length(Ks)
        for a = 1:length(alphas)
            for i = 1:length(resMVMD.imfs_MVMD)
                for j = 1:size(resMVMD.imfs_MVMD{i,k,a},1)
                    if i == 1 && j == 1
                        [~,f_mvmd] = pwelch(squeeze(resMVMD.imfs_MVMD{i,k,a}(j,:,:)),200,[],[],fs);
                    end
                    Px(i,j,:) = sum(pwelch(squeeze(resMVMD.imfs_MVMD{i,k,a}(j,:,:)),200,[],[],fs),2);
                end
            end
            for i = 1:size(Px,1)
                for j = 1:size(Px,2)
                    frequencies_MVMD_sep(i,j) = f_mvmd(Px(i,j,:) == max(Px(i,j,:)));
                end
            end
            Px = squeeze(sum(Px,1))';
            for i = 1:size(Px,2)
                frequencies_MVMD(i) = f_mvmd(Px(:,i) == max(Px(:,i)));
            end
            Pxs_MVMD = Px;
            [~,pos_maxs] = max(Pxs_MVMD);
            f_maxs = f_mvmd(pos_maxs);
            [~,inds_MVMD(k,a)] = min(abs(f_maxs-0.05)); %Find the component closest to 0.05 Hz
        end
    end
end
for k = 1:length(Ks)
    for a = 1:length(alphas)
        for i = 1:numSubjects
            if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_MVMD_PSA.mat']))
                COSDELPHI = phase_sync_analysis_HCP(resMVMD.imfs_MVMD{i,k,a},'MVMD',indx,...
                     inds_MVMD(k,a));
                COSDELPHI1_MVMD{i,k,a} = COSDELPHI;
            end
        end
    end
end

% if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\MSWD_HCP_PSA.mat']))
%     save(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\MSWD_HCP_PSA.mat']),"imfs_MSWD")
% end
% 
% if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_MSWD_PSA.mat']))
%     save(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_MSWD_PSA.mat']),'COSDELPHI1_MSWD','-v7.3')
% end


if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\MVMD_HCP_PSA.mat']))
    save(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\MVMD_HCP_PSA.mat']),"imfs_MVMD")
end

if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_MVMD_PSA.mat']))
    save(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_MVMD_PSA.mat']),'COSDELPHI1_MVMD','-v7.3')
end

if ~exist('COSDELPHI1_MSWD','var')
    load(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_MSWD_PSA.mat']))
end
if ~exist('COSDELPHI1_MVMD','var')
    load(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_MVMD_PSA.mat']))
end


DBI_MVMD = 2;
DBI_MSWD = 2;
for k = 1:length(Ks)
    for a = 1:length(alphas)
        for i = 1:length(COSDELPHI1_MSWD)
            if i == 1
                COSDELPHI_MSWD = COSDELPHI1_MSWD{i,k,a};
                COSDELPHI_MVMD = COSDELPHI1_MVMD{i,k,a};
            else
                COSDELPHI_MSWD = cat(1,COSDELPHI_MSWD,COSDELPHI1_MSWD{i,k,a});
                COSDELPHI_MVMD = cat(1,COSDELPHI_MVMD,COSDELPHI1_MVMD{i,k,a});
            end
        end
        if ~exist(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MVMD_PSA.mat'])) ||...
            ~exist(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MSWD_PSA.mat']))
            [~,C_init_MVMD,sumd_MVMD] = kmeans(squeeze(COSDELPHI_MVMD),DBI_MVMD,'MaxIter',150,'Start','sample','Replicates',10); % replicates default 200, MaxIter 150, replicates 10
            [mean_idx_MVMD_psa{k,a},mean_C_MVMD_psa{k,a}] = kmeans(squeeze(COSDELPHI_MVMD),DBI_MVMD,'MaxIter',500,'Start',C_init_MVMD); %1000 default
            
            [~,C_init_MSWD,sumd_MSWD] = kmeans(squeeze(COSDELPHI_MSWD),DBI_MSWD,'MaxIter',150,'Start','sample','Replicates',10); % replicates default 200, MaxIter 150, replicates 10
            [mean_idx_MSWD_psa{k,a},mean_C_MSWD_psa{k,a}] = kmeans(squeeze(COSDELPHI_MSWD),DBI_MSWD,'MaxIter',500,'Start',C_init_MSWD); %1000 default
            save(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MVMD_PSA.mat']),'mean_idx_MVMD_psa','mean_C_MVMD_psa');
            save(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MSWD_PSA.mat']),'mean_idx_MSWD_psa',"mean_C_MSWD_psa");
        else
            load(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MVMD_PSA.mat']));
            load(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MSWD_PSA.mat']));
        end

        load(fullfile(path_save,'MSWD_CL_paper_decomposed_HCP_rest1\MVMD_HCP_50_subjects.mat'))
        load(fullfile(path_save,'MSWD_CL_paper_decomposed_HCP_rest1\MSWD_HCP_50_subjects.mat'))

        r1 = corrcoef(mean_C_MVMD_psa{k,a}(1,:),mean_C_MVMD_all{1}(1,:));
        r2 = corrcoef(mean_C_MVMD_psa{k,a}(2,:),mean_C_MVMD_all{1}(2,:));
        r3 = corrcoef(mean_C_MVMD_psa{k,a}(1,:),mean_C_MVMD_all{1}(2,:));
        r4 = corrcoef(mean_C_MVMD_psa{k,a}(2,:),mean_C_MVMD_all{1}(1,:));
        if r1(1,2) + r2(1,2) > r3(1,2) + r4(1,2)
            corrs_MVMD(k,a,1) = r1(1,2);
            corrs_MVMD(k,a,2) = r2(1,2);
        else
            corrs_MVMD(k,a,1) = r4(1,2);
            corrs_MVMD(k,a,2) = r3(1,2);
        end
    
        r1 = corrcoef(mean_C_MSWD_psa{k,a}(1,:),mean_C_MSWD_all{1}(1,:));
        r2 = corrcoef(mean_C_MSWD_psa{k,a}(2,:),mean_C_MSWD_all{1}(2,:));
        r3 = corrcoef(mean_C_MSWD_psa{k,a}(1,:),mean_C_MSWD_all{1}(2,:));
        r4 = corrcoef(mean_C_MSWD_psa{k,a}(2,:),mean_C_MSWD_all{1}(1,:));
        if r1(1,2) + r2(1,2) > r3(1,2) + r4(1,2)
            corrs_MSWD(k,a,1) = r1(1,2);
            corrs_MSWD(k,a,2) = r2(1,2);
        else
            corrs_MSWD(k,a,1) = r4(1,2);
            corrs_MSWD(k,a,2) = r3(1,2);
        end 
    end
end


figure
surf(1:length(alphas),1:length(Ks),corrs_MVMD(:,:,1)');
ax = gca;
ax.FontSize = 12;
xlabel('K','FontSize',14)
xticks(Ks)
ylabel('alpha','FontSize',14)
yticks(1:length(alphas))
yticklabels({'500','2000','5000'})
zlim([minimum,maximum])
c=colorbar('Position',[1 0.001 0.001 0.01]); caxis([minimum+0.2 maximum-0.2])


figure
surf(1:length(P_corr_imps),1:length(compStds),corrs_MSWD(:,:,2)');
ax = gca;
ax.FontSize = 12;
xlabel('Corr_t_h','FontSize',14)
xticks(1:length(P_corr_imps))
ylabel('StD_t_h','FontSize',14)
yticks(1:length(compStds))
yticklabels({'0.002','0.01','0.1','0.2'})
xticklabels({'0.01','0.07','0.1'})
zlim([minimum,maximum])
c=colorbar('Position',[1 0.001 0.001 0.01]); caxis([minimum+0.2 maximum-0.2])


figure
surf(1:length(P_corr_imps),1:length(compStds),corrs_MSWD(:,:,1)');
ax = gca;
ax.FontSize = 12;
xlabel('Corr_t_h','FontSize',14)
xticks(1:length(P_corr_imps))
ylabel('StD_t_h','FontSize',14)
yticks(1:length(compStds))
yticklabels({'0.002','0.01','0.1','0.2'})
xticklabels({'0.01','0.07','0.1'})
zlim([minimum,maximum])
c=colorbar('Position',[1 0.001 0.001 0.01]); caxis([minimum+0.2 maximum-0.2])


figure
surf(1:length(P_corr_imps),1:length(compStds),corrs_MSWD(:,:,2)');
ax = gca;
ax.FontSize = 12;
xlabel('Corr_t_h','FontSize',14)
xticks(1:length(P_corr_imps))
ylabel('StD_t_h','FontSize',14)
yticks(1:length(compStds))
yticklabels({'0.002','0.01','0.1','0.2'})
xticklabels({'0.01','0.07','0.1'})
zlim([minimum,maximum])
c=colorbar('Position',[1 0.001 0.001 0.01]); caxis([minimum+0.2 maximum-0.2])

