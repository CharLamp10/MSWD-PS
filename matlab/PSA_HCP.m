% In this script, fMRI data from the HCP project are decomposed using MVMD
% and MSWD and CRP is calculated for dominant component

clear
clc

save_plots = true;
save_res = 1;
numICs = 100;
visit = 'rest1';
path_data = 'C:\Users\100063082\Desktop\MSWD_paper_files';
path_save = 'F:\MSWD_paper_new';

if contains(pwd,'100063082')
    path_tmaps = fullfile(path_data,'dyn_tmap_comps_centered.nii');
    path_hcp_average = fullfile(path_data,'HCP_PTN1200\groupICA\groupICA_3T_HCP1200_MSMAll_d100.ica\melodic_IC_sum.nii');
    path_data = fullfile(path_data,'HCP_PTN1200\NodeTimeseries_3T_HCP1200_MSMAll_ICAd100_ts2\node_timeseries\3T_HCP1200_MSMAll_d100_ts2');
end

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

numSubjects = 50;

if exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\MEMD_HCP_PSA.mat'])) &&...
        ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_MEMD_PSA.mat']))
    resMEMD = load(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\MEMD_HCP_PSA.mat']));
end
if exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\EWT_HCP_PSA.mat'])) &&...
        ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_EWT_PSA.mat']))
    resEWT = load(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\EWT_HCP_PSA.mat']));
end
if exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\MSWD_HCP_PSA.mat'])) &&...
        ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_MVMD_PSA.mat']))
    resMSWD = load(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\MSWD_HCP_PSA.mat']));
end
if exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\MVMD_HCP_PSA.mat'])) &&...
        ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_MSWD_PSA.mat']))
    resMVMD = load(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\MVMD_HCP_PSA.mat']));
end

load("names_98_subjects.mat");


inds_MSWD = [1]; %In MSWD we get always the first component
inds_MEMD = [1]; %We force the component with the maximum spectral peak within the [0.01,0.1] Hz band to be at the first index
inds_EWT = [1]; %We force the component with the maximum spectral peak within the [0.01,0.1] Hz band to be at the first index
for i = 1:numSubjects
    name = names{i};
    data = load(fullfile(path_data,name));
    data = data(1:1200,I);
    data = data./max(abs(data));
    if i == 1
        [~,f] = pwelch(data,size(data,1)/2,[],[],fs);
    end
    Px = mean(pwelch(data,size(data,1)/2,[],[],fs),2);
    [~,pos_max] = max(Px);
    f_EWT(i) = f(pos_max);
end

indx = nchoosek(1:numICs,2);

stop_vecs = {[0.075,0.75,0.075],[0.2,0.8,0.2],[0.3,0.3,0.3],[0.5,0.5,0.5]};
ndirs = [0.05,0.10,0.2];
P_ths = [0.1,0.3,0.5,0.8];
freq_ress = [0.08,0.16,0.25];
Ks = 1:2:7;
alphas = [500,2000,5000];
compStds = [0.002,0.01,0.1,0.2];
P_corr_imps = [0.01,0.07,0.1];

imfs_MEMD = cell(numSubjects,length(stop_vecs),length(ndirs));
imfs_EWT = cell(numSubjects,length(P_ths),length(freq_ress));
imfs_MSWD = cell(numSubjects,length(compStds),length(P_corr_imps));
imfs_MVMD = cell(numSubjects,length(Ks),length(alphas));

%% MEMD
for k = 1:length(stop_vecs)
    for a = 1:length(ndirs)
        for i = 1:numSubjects
            if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\MEMD_HCP_PSA.mat']))
                name = names{i};
                data = load(fullfile(path_data,name));
                data_MEMD = data(1:1200,I);
                data_MEMD = data_MEMD./max(abs(data_MEMD));
                stop_vec = stop_vecs{k};
                ndir = ndirs(a);
                stp_crit = 'stop';
                mode = 'na_snr';
                intensity_noise = 0.75; 
                n_channel_na = size(data_MEMD,2);  
                ndirr = ndir*n_channel_na;
                imfs = namemd(data_MEMD, ndirr, stp_crit, stop_vec, mode, intensity_noise, n_channel_na);
                for j = 1:length(imfs)
                    imf(:,:,j) = imfs{j};
                end
                for j = 1:size(imf,1)
                    p(:,j) = mean(pwelch(squeeze(imf(j,:,:)),size(imf,2),[],[],fs),2);
                end
                f = pwelch(squeeze(imf(1,:,1)),size(imf,2),[],[],fs);
                [maxs,pos_maxs] = max(p);
                clear p
                maxi = 0;
                for j = 1:length(maxs)
                    if maxs(j)> maxi && f(pos_maxs(j)) >= 0.01 && f(pos_maxs(j)) <= 0.1
                        maxi = maxs(j);
                    end
                end
                pos_max = find(maxs == maxi);
                temp = imf(1,:,:);
                imf(1,:,:) = imf(pos_max,:,:);
                imf(pos_max,:,:) = temp;
                imfs_MEMD{i,k,a} = imf;
                disp(['ndir:', num2str(ndirs(a)), ' stop_vec:', num2str(stop_vecs{k})])
            elseif ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_MEMD_PSA.mat']))
                imf = resMEMD.imfs_MEMD{i,k,a};
            end
            if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_MEMD_PSA.mat']))
                COSDELPHI = phase_sync_analysis_HCP(imf,'MVMD',indx,inds_MSWD);
                COSDELPHI1_MEMD{i,k,a} = COSDELPHI;
            end
            clear imf
        end
    end
end

%% EWT
for k = 1:length(P_ths)
    for a = 1:length(freq_ress)
        for i = 1:numSubjects
            if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\EWT_HCP_PSA.mat']))
                name = names{i};
                data = load(fullfile(path_data,name));
                data_EWT = data(1:1200,I);
                data_EWT = data_EWT./max(abs(data_EWT));
                P_th = P_ths(k);
                freq_res = freq_ress(a);
                imf = zeros(20,size(data_EWT,1),size(data_EWT,2));
                for j = 1:size(data_EWT,2)
                    temp_imf = ewt(data_EWT(:,j),"PeakThresholdPercent",P_th,"FrequencyResolution",freq_res)';
                    mfreqs = zeros(1,size(temp_imf,1));
                    for n = 1:size(temp_imf,1)
                        mfreqs(n) = meanfreq(temp_imf(n,:),fs);
                    end
                    [~,pos] = min(abs(mfreqs - f_EWT(i)));
                    temp = temp_imf(pos,:);
                    temp_imf(pos,:) = temp_imf(1,:);
                    temp_imf(1,:) = temp;
                    imf(1:size(temp_imf,1),:,j) = temp_imf;
                    clear temp_imf
                end
                imf(all(all(imf == 0,3),2),:,:) = [];
                imfs_EWT{i,k,a} = imf;
                disp(['P_th:', num2str(P_ths(k)), ' freq_res:', num2str(freq_ress(a))])
            elseif ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_EWT_PSA.mat']))
                imf = resEWT.imfs_EWT{i,k,a};
            end
            if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_EWT_PSA.mat']))
                COSDELPHI = phase_sync_analysis_HCP(imf,'MVMD',indx,inds_EWT);
                COSDELPHI1_EWT{i,k,a} = COSDELPHI;
            end
            clear imf
        end
    end
end

%% MVMD
for k = 1:length(Ks)
    for a = 1:length(alphas)
        for i = 1:numSubjects
            if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\MVMD_HCP_PSA.mat']))
                name = names{i};
                data = load(fullfile(path_data,name));
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
                        [~,f_mvmd] = pwelch(squeeze(resMVMD.imfs_MVMD{i,k,a}(j,:,:)),1200,[],[],fs);
                    end
                    Px(i,j,:) = sum(pwelch(squeeze(resMVMD.imfs_MVMD{i,k,a}(j,:,:)),1200,[],[],fs),2);
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
            clear Px
        end
    end

    for k = 1:length(Ks)
        for a = 1:length(alphas)
            for i = 1:numSubjects
                if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_MVMD_PSA.mat']))
                    COSDELPHI = phase_sync_analysis_HCP(resMVMD.imfs_MVMD{i,k,a},'MVMD',indx,inds_MVMD(k,a));
                    COSDELPHI1_MVMD{i,k,a} = COSDELPHI;
                end
            end
        end
    end
end

%% MSWD
for k = 1:length(compStds)
    for a = 1:length(P_corr_imps)
        for i = 1:numSubjects
            if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\MSWD_HCP_PSA.mat']))
                name = names{i};
                data = load(fullfile(path_data,name));
                data_MSWD = data(1:1200,I);
                data_MSWD = data_MSWD./max(abs(data_MSWD));
                p_value = 1e-5;
                P_corr = 1;
                P_corr_imp = P_corr_imps(a);
                wind = [];
                compStd = compStds(k);
                param_struct  = struct('P_corr', P_corr, ...
                    'P_corr_imp',   P_corr_imp,...
                    'StD_th',       compStd, ...
                    'Welch_window', wind, ...
                    'p_value',      1e-5);
                imf = MSWD(data_MSWD, param_struct);
                imfs_MSWD{i,k,a} = imf;
                disp(['P_corr_imp:', num2str(P_corr_imps(a)), ' compStd:', num2str(compStds(k))])
            elseif ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_MSWD_PSA.mat']))
                imf = resMSWD.imfs_MSWD{i,k,a};
            end
            if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_MSWD_PSA.mat']))
                COSDELPHI = phase_sync_analysis_HCP(imf,'MSWD',indx,inds_MSWD);
                COSDELPHI1_MSWD{i,k,a} = COSDELPHI;
            end
        end
    end
end


if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\MEMD_HCP_PSA.mat']))
    save(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\MEMD_HCP_PSA.mat']),"imfs_MEMD",'-v7.3')
end

if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\EWT_HCP_PSA.mat']))
    save(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\EWT_HCP_PSA.mat']),"imfs_EWT")
end

if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\MVMD_HCP_PSA.mat']))
    save(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\MVMD_HCP_PSA.mat']),"imfs_MVMD")
end

if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\MSWD_HCP_PSA.mat']))
    save(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\MSWD_HCP_PSA.mat']),"imfs_MSWD")
end

if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_MEMD_PSA.mat']))
    save(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_MEMD_PSA.mat']),'COSDELPHI1_MEMD','-v7.3')
end

if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_EWT_PSA.mat']))
    save(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_EWT_PSA.mat']),'COSDELPHI1_EWT','-v7.3')
end

if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_MSWD_PSA.mat']))
    save(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_MSWD_PSA.mat']),'COSDELPHI1_MSWD','-v7.3')
end

if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_MVMD_PSA.mat']))
    save(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_MVMD_PSA.mat']),'COSDELPHI1_MVMD','-v7.3')
end

if ~exist('COSDELPHI1_MEMD','var') && ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MEMD_PSA.mat']))
    load(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_MEMD_PSA.mat']))
end
if ~exist('COSDELPHI1_EWT','var') && ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_EWT_PSA.mat']))
    load(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_EWT_PSA.mat']))
end
if ~exist('COSDELPHI1_MVMD','var') && ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MVMD_PSA.mat']))
    load(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_MVMD_PSA.mat']))
end
if ~exist('COSDELPHI1_MSWD','var') && ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MSWD_PSA.mat']))
    load(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_MSWD_PSA.mat']))
end

DBI_MEMD = 2;
DBI_EWT = 2;
DBI_MVMD = 2;
DBI_MSWD = 2;
 if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MEMD_PSA.mat'])) ||...
      ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_EWT_PSA.mat'])) ||...
       ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MVMD_PSA.mat'])) ||...
        ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MSWD_PSA.mat']))
    for k = 1:length(Ks)
        for a = 1:length(alphas)
            for i = 1:length(COSDELPHI1_MSWD)
                if i == 1
                    COSDELPHI_MEMD = COSDELPHI1_MEMD{i,k,a};
                    COSDELPHI_EWT = COSDELPHI1_EWT{i,k,a};
                    COSDELPHI_MSWD = COSDELPHI1_MSWD{i,k,a};
                    COSDELPHI_MVMD = COSDELPHI1_MVMD{i,k,a};
                else
                    COSDELPHI_MEMD = cat(1,COSDELPHI_MEMD,COSDELPHI1_MEMD{i,k,a});
                    COSDELPHI_EWT = cat(1,COSDELPHI_EWT,COSDELPHI1_EWT{i,k,a});
                    COSDELPHI_MSWD = cat(1,COSDELPHI_MSWD,COSDELPHI1_MSWD{i,k,a});
                    COSDELPHI_MVMD = cat(1,COSDELPHI_MVMD,COSDELPHI1_MVMD{i,k,a});
                end
            end
            if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MEMD_PSA.mat']))
                [~,C_init_MEMD,sumd_MEMD] = kmeans(squeeze(COSDELPHI_MEMD),DBI_MEMD,'MaxIter',150,'Start','sample','Replicates',10); % replicates default 200, MaxIter 150, replicates 10
                [mean_idx_MEMD_psa{k,a},mean_C_MEMD_psa{k,a}] = kmeans(squeeze(COSDELPHI_MEMD),DBI_MEMD,'MaxIter',500,'Start',C_init_MEMD); %1000 default
            end
            if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_EWT_PSA.mat']))
                [~,C_init_EWT,sumd_EWT] = kmeans(squeeze(COSDELPHI_EWT),DBI_EWT,'MaxIter',150,'Start','sample','Replicates',10); % replicates default 200, MaxIter 150, replicates 10
                [mean_idx_EWT_psa{k,a},mean_C_EWT_psa{k,a}] = kmeans(squeeze(COSDELPHI_EWT),DBI_EWT,'MaxIter',500,'Start',C_init_EWT); %1000 default
            end
            if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MVMD_PSA.mat']))
                [~,C_init_MVMD,sumd_MVMD] = kmeans(squeeze(COSDELPHI_MVMD),DBI_MVMD,'MaxIter',150,'Start','sample','Replicates',10); % replicates default 200, MaxIter 150, replicates 10
                [mean_idx_MVMD_psa{k,a},mean_C_MVMD_psa{k,a}] = kmeans(squeeze(COSDELPHI_MVMD),DBI_MVMD,'MaxIter',500,'Start',C_init_MVMD); %1000 default
            end
            if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MSWD_PSA.mat']))
                [~,C_init_MSWD,sumd_MSWD] = kmeans(squeeze(COSDELPHI_MSWD),DBI_MSWD,'MaxIter',150,'Start','sample','Replicates',10); % replicates default 200, MaxIter 150, replicates 10
                [mean_idx_MSWD_psa{k,a},mean_C_MSWD_psa{k,a}] = kmeans(squeeze(COSDELPHI_MSWD),DBI_MSWD,'MaxIter',500,'Start',C_init_MSWD); %1000 default
            end
            disp(['k = ', num2str(k), ', a = ', num2str(a)])
        end
    end
    if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MEMD_PSA.mat']))
        save(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MEMD_PSA.mat']),'mean_C_MEMD_psa','mean_idx_MEMD_psa');
    end
    if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_EWT_PSA.mat']))
        save(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_EWT_PSA.mat']),'mean_C_EWT_psa','mean_idx_EWT_psa');
    end
    if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MVMD_PSA.mat']))
        save(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MVMD_PSA.mat']),'mean_C_MVMD_psa','mean_idx_MVMD_psa');
    end
    if ~exist(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MSWD_PSA.mat']))
        save(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MSWD_PSA.mat']),'mean_C_MSWD_psa','mean_idx_MSWD_psa');
    end
 else
    if ~exist("mean_idx_MEMD_psa")
        load(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MEMD_PSA.mat']));
    end
    if ~exist("mean_idx_EWT_psa")
        load(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_EWT_PSA.mat']));
    end
    if ~exist("mean_idx_MVMD_psa")
        load(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MVMD_PSA.mat']));
    end
    if ~exist("mean_idx_MSWD_psa")
        load(fullfile(path_save,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MSWD_PSA.mat']));
    end
end

load(fullfile(path_save,'MSWD_CL_paper_decomposed_HCP_rest1\mean_C_MEMD_all.mat'))
load(fullfile(path_save,'MSWD_CL_paper_decomposed_HCP_rest1\mean_C_EWT_all.mat'))
load(fullfile(path_save,'MSWD_CL_paper_decomposed_HCP_rest1\mean_C_MVMD_all.mat'))
load(fullfile(path_save,'MSWD_CL_paper_decomposed_HCP_rest1\mean_C_MSWD_all.mat'))

for k = 1:length(Ks)
    for a = 1:length(alphas)
        r1 = rmse(mean_C_MEMD_psa{k,a}(1,:),mean_C_MEMD_all{1}(1,:));
        r2 = rmse(mean_C_MEMD_psa{k,a}(2,:),mean_C_MEMD_all{1}(2,:));
        r3 = rmse(mean_C_MEMD_psa{k,a}(1,:),mean_C_MEMD_all{1}(2,:));
        r4 = rmse(mean_C_MEMD_psa{k,a}(2,:),mean_C_MEMD_all{1}(1,:));
        if r1 + r2 < r3 + r4
            corrs_MEMD(k,a,1) = r1;
            corrs_MEMD(k,a,2) = r2;
        else
            corrs_MEMD(k,a,1) = r4;
            corrs_MEMD(k,a,2) = r3;
        end

        r1 = rmse(mean_C_EWT_psa{k,a}(1,:),mean_C_EWT_all{1}(1,:));
        r2 = rmse(mean_C_EWT_psa{k,a}(2,:),mean_C_EWT_all{1}(2,:));
        r3 = rmse(mean_C_EWT_psa{k,a}(1,:),mean_C_EWT_all{1}(2,:));
        r4 = rmse(mean_C_EWT_psa{k,a}(2,:),mean_C_EWT_all{1}(1,:));
        if r1 + r2 < r3 + r4
            corrs_EWT(k,a,1) = r1;
            corrs_EWT(k,a,2) = r2;
        else
            corrs_EWT(k,a,1) = r4;
            corrs_EWT(k,a,2) = r3;
        end

        r1 = rmse(mean_C_MVMD_psa{k,a}(1,:),mean_C_MVMD_all{1}(1,:));
        r2 = rmse(mean_C_MVMD_psa{k,a}(2,:),mean_C_MVMD_all{1}(2,:));
        r3 = rmse(mean_C_MVMD_psa{k,a}(1,:),mean_C_MVMD_all{1}(2,:));
        r4 = rmse(mean_C_MVMD_psa{k,a}(2,:),mean_C_MVMD_all{1}(1,:));
        if r1 + r2 < r3 + r4
            corrs_MVMD(k,a,1) = r1;
            corrs_MVMD(k,a,2) = r2;
        else
            corrs_MVMD(k,a,1) = r4;
            corrs_MVMD(k,a,2) = r3;
        end
    
        r1 = rmse(mean_C_MSWD_psa{k,a}(1,:),mean_C_MSWD_all{1}(1,:));
        r2 = rmse(mean_C_MSWD_psa{k,a}(2,:),mean_C_MSWD_all{1}(2,:));
        r3 = rmse(mean_C_MSWD_psa{k,a}(1,:),mean_C_MSWD_all{1}(2,:));
        r4 = rmse(mean_C_MSWD_psa{k,a}(2,:),mean_C_MSWD_all{1}(1,:));
        if r1 + r2 < r3 + r4
            corrs_MSWD(k,a,1) = r1;
            corrs_MSWD(k,a,2) = r2;
        else
            corrs_MSWD(k,a,1) = r4;
            corrs_MSWD(k,a,2) = r3;
        end 
    end
end

%% MVMD-MSWD
minimum = min([min(corrs_MVMD(:,:,1),[],"all"),min(corrs_MSWD(:,:,1),[],"all")]);
maximum = max([max(corrs_MVMD(:,:,1),[],"all"),max(corrs_MSWD(:,:,1),[],"all")]);

figure
surf(1:length(alphas),Ks,corrs_MVMD(:,:,1));
ax = gca;
ax.FontSize = 12;
ylabel('K','FontSize',14)
yticks(Ks)
xlabel('alpha','FontSize',14)
xticks(1:length(alphas))
xticklabels({'500','2000','5000'})
zlim([minimum,maximum])
c=colorbar('Position',[1 0.001 0.001 0.01]); caxis([minimum maximum])
if save_plots
    exportgraphics(gcf,fullfile('plots_PSA_new','MVMD_state_1_HCP_RMSE.png'),'Resolution',1000)
end

figure
surf(1:length(compStds),1:length(P_corr_imps),corrs_MSWD(:,:,1)');
ax = gca;
ax.FontSize = 12;
ylabel('Corr_t_h','FontSize',14)
yticks(1:length(P_corr_imps))
xlabel('StD_t_h','FontSize',14)
xticks(1:length(compStds))
xticklabels({'0.002','0.01','0.1','0.2'})
yticklabels({'0.01','0.07','0.1'})
zlim([minimum,maximum])
c=colorbar('Position',[1 0.001 0.001 0.01]); caxis([minimum maximum])
if save_plots
    exportgraphics(gcf,fullfile('plots_PSA_new','MSWD_state_1_HCP_RMSE.png'),'Resolution',1000)
end


minimum = min([min(corrs_MVMD(:,:,2),[],"all"),min(corrs_MSWD(:,:,2),[],"all")]);
maximum = max([max(corrs_MVMD(:,:,2),[],"all"),max(corrs_MSWD(:,:,2),[],"all")]);

figure
surf(1:length(alphas),Ks,corrs_MVMD(:,:,2));
ax = gca;
ax.FontSize = 12;
ylabel('K','FontSize',14)
yticks(Ks)
xlabel('alpha','FontSize',14)
xticks(1:length(alphas))
xticklabels({'500','2000','5000'})
zlim([minimum,maximum])
c=colorbar('Position',[1 0.001 0.001 0.01]); caxis([minimum maximum])
if save_plots
    exportgraphics(gcf,fullfile('plots_PSA_new','MVMD_state_2_HCP_RMSE.png'),'Resolution',1000)
end


figure
surf(1:length(compStds),1:length(P_corr_imps),corrs_MSWD(:,:,2)');
ax = gca;
ax.FontSize = 12;
ylabel('Corr_t_h','FontSize',14)
yticks(1:length(P_corr_imps))
xlabel('StD_t_h','FontSize',14)
xticks(1:length(compStds))
xticklabels({'0.002','0.01','0.1','0.2'})
yticklabels({'0.01','0.07','0.1'})
zlim([minimum,maximum])
c=colorbar('Position',[1 0.001 0.001 0.01]); caxis([minimum maximum])
if save_plots
    exportgraphics(gcf,fullfile('plots_PSA_new','MSWD_state_2_HCP_RMSE.png'),'Resolution',1000)
end



%% MEMD-EWT-MVMD-MSWD
minimum = min([min(corrs_MEMD(:,:,1),[],"all"),min(corrs_EWT(:,:,1),[],"all"),min(corrs_MVMD(:,:,1),[],"all"),min(corrs_MSWD(:,:,1),[],"all")]);
maximum = max([max(corrs_MEMD(:,:,1),[],"all"),max(corrs_EWT(:,:,1),[],"all"),max(corrs_MVMD(:,:,1),[],"all"),max(corrs_MSWD(:,:,1),[],"all")]);

figure
surf(1:length(ndirs),1:length(stop_vecs),corrs_MEMD(:,:,1));
ax = gca;
ax.FontSize = 12;
xlabel('$p$','Interpreter','latex','FontSize',18)
yticks(1:length(stop_vecs))
yticklabels({'[0.075,0.75,0.075]','[0.2,0.9,0.2]','[0.3,0.3,0.3]','[0.5,0.5,0.5]'})
ylabel('[$sd_1$,$sd_2$,$tol$]','Interpreter','latex','FontSize',18)
xticks(1:length(ndirs))
xticklabels({'5','10','20'})
zlim([minimum,maximum])
c=colorbar('Position',[1 0.001 0.001 0.01]); caxis([minimum maximum])
if save_plots
    exportgraphics(gcf,fullfile('plots_PSA_new','MEMD_state_1_HCP_RMSE.png'),'Resolution',1000)
end

figure
surf(1:length(freq_ress),1:length(P_ths),corrs_EWT(:,:,1));
ax = gca;
ax.FontSize = 12;
ylabel('P_t_h','FontSize',14)
yticks(1:length(P_ths))
xlabel('freqRes','FontSize',14)
xticks(1:length(freq_ress))
xticklabels({'0.08','0.16','0.25'}) 
yticklabels({'0.1','0.3','0.5','0.8'})
zlim([minimum,maximum])
c=colorbar('Position',[1 0.001 0.001 0.01]); caxis([minimum maximum])
if save_plots
    exportgraphics(gcf,fullfile('plots_PSA_new','EWT_state_1_HCP_RMSE.png'),'Resolution',1000)
end

figure
surf(1:length(alphas),Ks,corrs_MVMD(:,:,1));
ax = gca;
ax.FontSize = 12;
ylabel('K','FontSize',14)
yticks(Ks)
xlabel('alpha','FontSize',14)
xticks(1:length(alphas))
xticklabels({'500','2000','5000'})
zlim([minimum,maximum])
c=colorbar('Position',[1 0.001 0.001 0.01]); caxis([minimum maximum])
if save_plots
    exportgraphics(gcf,fullfile('plots_PSA_new','MVMD_all_state_1_HCP_RMSE.png'),'Resolution',1000)
end

figure
surf(1:length(compStds),1:length(P_corr_imps),corrs_MSWD(:,:,1)');
ax = gca;
ax.FontSize = 12;
ylabel('Corr_t_h','FontSize',14)
yticks(1:length(P_corr_imps))
xlabel('StD_t_h','FontSize',14)
xticks(1:length(compStds))
xticklabels({'0.002','0.01','0.1','0.2'})
yticklabels({'0.01','0.07','0.1'})
zlim([minimum,maximum])
c=colorbar('Position',[1 0.001 0.001 0.01]); caxis([minimum maximum])
if save_plots
    exportgraphics(gcf,fullfile('plots_PSA_new','MSWD_all_state_1_HCP_RMSE.png'),'Resolution',1000)
end


minimum = min([min(corrs_MEMD(:,:,2),[],"all"),min(corrs_EWT(:,:,2),[],"all"),min(corrs_MVMD(:,:,2),[],"all"),min(corrs_MSWD(:,:,2),[],"all")]);
maximum = max([max(corrs_MEMD(:,:,2),[],"all"),max(corrs_EWT(:,:,2),[],"all"),max(corrs_MVMD(:,:,2),[],"all"),max(corrs_MSWD(:,:,2),[],"all")]);

figure
surf(1:length(ndirs),1:length(stop_vecs),corrs_MEMD(:,:,2));
ax = gca;
ax.FontSize = 12;
xlabel('$p$','Interpreter','latex','FontSize',18)
yticks(1:length(stop_vecs))
yticklabels({'[0.075,0.75,0.075]','[0.2,0.9,0.2]','[0.3,0.3,0.3]','[0.5,0.5,0.5]'})
ylabel('[$sd_1$,$sd_2$,$tol$]','Interpreter','latex','FontSize',18)
xticks(1:length(ndirs))
xticklabels({'5','10','20'})
zlim([minimum,maximum])
c=colorbar('Position',[1 0.001 0.001 0.01]); caxis([minimum maximum])
if save_plots
    exportgraphics(gcf,fullfile('plots_PSA_new','MEMD_state_2_HCP_RMSE.png'),'Resolution',1000)
end


figure
surf(1:length(freq_ress),1:length(P_ths),corrs_EWT(:,:,2));
ax = gca;
ax.FontSize = 12;
ylabel('P_t_h','FontSize',14)
yticks(1:length(P_ths))
xlabel('freqRes','FontSize',14)
xticks(1:length(freq_ress))
xticklabels({'0.08','0.16','0.25'}) 
yticklabels({'0.1','0.3','0.5','0.8'})
zlim([minimum,maximum])
c=colorbar('Position',[1 0.001 0.001 0.01]); caxis([minimum maximum])
if save_plots
    exportgraphics(gcf,fullfile('plots_PSA_new','EWT_state_2_HCP_RMSE.png'),'Resolution',1000)
end



figure
surf(1:length(alphas),Ks,corrs_MVMD(:,:,2));
ax = gca;
ax.FontSize = 12;
ylabel('K','FontSize',14)
yticks(Ks)
xlabel('alpha','FontSize',14)
xticks(1:length(alphas))
xticklabels({'500','2000','5000'})
zlim([minimum-0.01,maximum])
c=colorbar('Position',[1 0.001 0.001 0.01]); caxis([minimum maximum])
if save_plots
    exportgraphics(gcf,fullfile('plots_PSA_new','MVMD_all_state_2_HCP_RMSE.png'),'Resolution',1000)
end


figure
surf(1:length(compStds),1:length(P_corr_imps),corrs_MSWD(:,:,2)');
ax = gca;
ax.FontSize = 12;
ylabel('Corr_t_h','FontSize',14)
yticks(1:length(P_corr_imps))
xlabel('StD_t_h','FontSize',14)
xticks(1:length(compStds))
xticklabels({'0.002','0.01','0.1','0.2'})
yticklabels({'0.01','0.07','0.1'})
zlim([minimum,maximum])
c=colorbar('Position',[1 0.001 0.001 0.01]); caxis([minimum maximum])
if save_plots
    exportgraphics(gcf,fullfile('plots_PSA_new','MSWD_all_state_2_HCP_RMSE.png'),'Resolution',1000)
end

