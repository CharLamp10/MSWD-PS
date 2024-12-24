% In this script, fMRI data from the HCP project are decomposed using MVMD
% and MSWD and CRP is calculated for dominant component

clear
clc

save_res = 1;
numICs = 100;
visit = 'rest1';
path_data = 'C:\Users\100063082\Desktop\MSWD_paper_files';
path_save = 'E:\MSWD_paper_new';

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

if exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\MSWD_HCP_50_subjects.mat']))
    resMSWD = load(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\MSWD_HCP_50_subjects.mat']));
end
if exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\MVMD_HCP_50_subjects.mat']))
    resMVMD = load(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\MVMD_HCP_50_subjects.mat']));
end

load("names_98_subjects.mat");


inds_MSWD = [1]; %In MSWD we get always the first component

indx = nchoosek(1:numICs,2);

for i = 1:numSubjects
    name = names{i};
    data = load(fullfile(path_data,name));
    %% MSWD
    if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\MSWD_HCP_98_subjects.mat']))
        data_MSWD = data(1:1200,I);
        data_MSWD = data_MSWD./max(abs(data_MSWD));
        p_value = 1e-5;
        P_corr = 1;
        P_corr_imp = 0.02;
        wind = [];
        compStd = 0.05;
        param_struct  = struct('P_corr', P_corr, ...
            'P_corr_imp',   P_corr_imp,...
            'StD_th',       compStd, ...
            'Welch_window', wind, ...
            'p_value',      1e-5);
        imf = MSWD(data_MSWD, param_struct);
        imfs_MSWD{i} = imf;
        disp(i)
    elseif ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_MSWD_50_subjects.mat']))
        imf = resMSWD.imfs_MSWD{i};
    end
    if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_MSWD_50_subjects.mat'])) && i <= 50
        COSDELPHI = phase_sync_analysis_HCP(imf,'MSWD',indx,...
             inds_MSWD);
        COSDELPHI1_MSWD{i} = COSDELPHI;
    end
end
    
for i = 1:numSubjects
    name = names{i};
    data = load(fullfile(path_data,name));
    %% MVMD
    if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\MVMD_HCP_50_subjects.mat']))
        data_MVMD = data(1:1200,I);
        data_MVMD = data_MVMD./max(abs(data_MVMD));
        tau = 0; DC = 1; init = 0; tol = 1e-9;
        K = 10; alpha = 1000;
        imf = MVMD_new(data_MVMD,alpha,tau,K,DC,init,tol);
        imfs_MVMD{i} = imf;
        disp(i)
    end
end
if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\MVMD_HCP_50_subjects.mat']))
    resMVMD.imfs_MVMD = imfs_MVMD;
end
if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_MVMD_50_subjects.mat']))
    for i = 1:length(resMVMD.imfs_MVMD)
        for j = 1:size(resMVMD.imfs_MVMD{i},1)
            if i == 1 && j == 1
                [~,f_mvmd] = pwelch(squeeze(resMVMD.imfs_MVMD{i}(j,:,:)),1200,[],[],fs);
            end
            Px(i,j,:) = sum(pwelch(squeeze(resMVMD.imfs_MVMD{i}(j,:,:)),1200,[],[],fs),2);
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
    [~,inds_MVMD] = min(abs(f_maxs-0.05)); %Find the component closest to 0.05 Hz
end
for i = 1:numSubjects
    if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_MVMD_50_subjects.mat']))
        COSDELPHI = phase_sync_analysis_HCP(resMVMD.imfs_MVMD{i},'MVMD',indx,...
             inds_MVMD);
        COSDELPHI1_MVMD{i} = COSDELPHI;
    end
end

if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\MSWD_HCP_98_subjects.mat']))
    save(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\MSWD_HCP_98_subjects.mat']),"imfs_MSWD")
end

if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_MSWD_50_subjects.mat']))
    save(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_MSWD_50_subjects.mat']),'COSDELPHI1_MSWD','-v7.3')
end


if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\MVMD_HCP_50_subjects.mat']))
    save(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\MVMD_HCP_50_subjects.mat']),"imfs_MVMD")
end

if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_MVMD_50_subjects.mat']))
    save(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_MVMD_50_subjects.mat']),'COSDELPHI1_MVMD','-v7.3')
end



